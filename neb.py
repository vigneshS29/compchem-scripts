#!/usr/bin/env python3
import os
import sys
import warnings
from copy import deepcopy
from collections import defaultdict, Counter

import numpy as np

import torch
import torch._dynamo
torch._dynamo.config.suppress_errors = True  # ok to keep; beware it can hide real errors

from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize import LBFGS
from ase.mep import NEB
from ase.data import covalent_radii, atomic_masses, chemical_symbols

# --- FAIR-Chem ---
from fairchem.core import pretrained_mlip, FAIRChemCalculator

import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment
import argparse

# =======================
# Helpers: FAIR-Chem
# =======================
def make_fairchem_predictor(model_name, device):
    return pretrained_mlip.get_predict_unit(model_name, device=device)

def attach_fairchem_calc(atoms, predictor):
    atoms.calc = FAIRChemCalculator(predictor, task_name="omol")

# =======================
# Geometry / graph helpers
# =======================
def _pairwise_distances(positions, cell=None, pbc=None):
    N = len(positions)
    if cell is None or pbc is None or not np.any(pbc):
        diffs = positions[:, None, :] - positions[None, :, :]
        return np.linalg.norm(diffs, axis=2)
    frac = cell.cartesian_to_fractional(positions)
    dfrac = frac[:, None, :] - frac[None, :, :]
    dfrac -= np.round(dfrac)
    dcart = cell.fractional_to_cartesian(dfrac.reshape(-1, 3)).reshape(N, N, 3)
    return np.linalg.norm(dcart, axis=2)

def _adjacency_by_covalent(atoms, cov_scale):
    Z = atoms.numbers
    R = covalent_radii[Z]
    D = _pairwise_distances(atoms.get_positions(), atoms.cell, atoms.pbc)
    cutoff = cov_scale * (R[:, None] + R[None, :])
    A = (D < cutoff) & (D > 1e-8)
    A = np.logical_and(A, A.T)
    np.fill_diagonal(A, False)
    return A, D

# =======================
# Depth-2 context weighting
# =======================
def _pair_sym(z1, z2):
    return tuple(sorted((chemical_symbols[int(z1)], chemical_symbols[int(z2)])))

def _node_base(z, base_mode):
    return float(int(z)) if base_mode == "Z" else float(atomic_masses[int(z)])

def _depth2_node_score(i, A, Z, base_mode, exclude=None, d1c=1.0, d2c=0.5):
    base_i = _node_base(Z[i], base_mode)
    n1 = np.where(A[i])[0].tolist()
    if exclude is not None and exclude in n1:
        n1.remove(exclude)
    s1 = sum(_node_base(Z[k], base_mode) for k in n1)
    n2_set = set()
    for k in n1:
        for l in np.where(A[k])[0]:
            if l == i or l == exclude or l in n1:
                continue
            n2_set.add(int(l))
    s2 = sum(_node_base(Z[l], base_mode) for l in n2_set)
    return base_i + d1c * s1 + d2c * s2

def _pair_context_weight(i, j, A, Z, base_mode, pair_mult_sym, bond_weight_exp):
    si = _depth2_node_score(i, A, Z, base_mode, exclude=j)
    sj = _depth2_node_score(j, A, Z, base_mode, exclude=i)
    w = (si * sj) ** 0.5
    if bond_weight_exp != 1.0:
        w = w ** bond_weight_exp
    mult = pair_mult_sym.get(_pair_sym(Z[i], Z[j]), 1.0)
    return w * mult

def weight_matrix_depth2(atoms, A, base_mode="mass", bond_weight_exp=1.0, pair_mult_sym=None):
    if pair_mult_sym is None:
        pair_mult_sym = {}
    Z = atoms.numbers.astype(int)
    n = len(Z)
    W = np.zeros((n, n), float)
    for i in range(n):
        for j in np.where(A[i])[0]:
            if j <= i:
                continue
            w = _pair_context_weight(i, j, A, Z, base_mode, pair_mult_sym, bond_weight_exp)
            W[i, j] = W[j, i] = w
    return W

# =======================
# Mapping costs (Hungarian seed)
# =======================
def _atom_cost_endpoint(i, j, A_m, A_t, Dm, Dt, Wm, Wt, Zm, Zt, k_dist=4):
    deg_m = float(np.sum(Wm[i, A_m[i]])) if np.any(A_m[i]) else 0.0
    deg_t = float(np.sum(Wt[j, A_t[j]])) if np.any(A_t[j]) else 0.0
    c_deg = abs(deg_m - deg_t)

    neigh_m = Counter(int(Zm[k]) for k in np.where(A_m[i])[0])
    neigh_t = Counter(int(Zt[k]) for k in np.where(A_t[j])[0])
    allZ = set(neigh_m) | set(neigh_t)
    c_neigh = sum(abs(neigh_m.get(z, 0) - neigh_t.get(z, 0)) for z in allZ)

    di = np.sort(Dm[i, A_m[i]])[:k_dist]
    dj = np.sort(Dt[j, A_t[j]])[:k_dist]
    if di.size < k_dist: di = np.pad(di, (0, k_dist - di.size))
    if dj.size < k_dist: dj = np.pad(dj, (0, k_dist - dj.size))
    c_dist = np.linalg.norm(di - dj)

    return 3.0 * c_deg + 1.0 * c_neigh + 1.0 * c_dist

# =======================
# Kabsch + RMSD
# =======================
def kabsch_align(mobile_atoms, target_atoms, use_masses=True):
    mob = mobile_atoms.copy()
    X = mob.get_positions()
    Y = target_atoms.get_positions()
    w = mob.get_masses() if use_masses else np.ones(len(mob))
    w = w / w.sum()
    Xc = X - np.average(X, axis=0, weights=w)
    Yc = Y - np.average(Y, axis=0, weights=w)
    H = (Xc * w[:, None]).T @ Yc
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    X_aligned = Xc @ R + np.average(Y, axis=0, weights=w)
    mob.set_positions(X_aligned)
    mob.set_cell(target_atoms.cell); mob.set_pbc(target_atoms.pbc)
    return mob

def rmsd(a, b):
    dx = a.get_positions() - b.get_positions()
    return float(np.sqrt((dx * dx).sum() / len(a)))

# =======================
# Breaks-only objective (vs initial)
# =======================
def eval_breaks_only(perm, A_init, A_t, W_init):
    breaks_w = 0.0
    keep_w   = 0.0
    br = 0
    keep = 0
    n = len(perm)
    for i in range(n):
        pi = perm[i]
        for j in range(i+1, n):
            pj = perm[j]
            init_bond = A_init[pi, pj]
            img_bond  = A_t[i, j]
            if init_bond:
                if img_bond:
                    keep_w += W_init[pi, pj]; keep += 1
                else:
                    breaks_w += W_init[pi, pj]; br   += 1
    return keep_w, breaks_w, keep, br

def total_breaks_across_images(perm, A_init, W_init, A_list):
    tot_keep_w = 0.0
    tot_breaks_w = 0.0
    tot_keep = 0
    tot_breaks = 0
    for A_t in A_list:
        keep_w, breaks_w, keep, br = eval_breaks_only(perm, A_init, A_t, W_init)
        tot_keep_w   += keep_w
        tot_breaks_w += breaks_w
        tot_keep     += keep
        tot_breaks   += br
    return tot_keep_w, tot_breaks_w, tot_keep, tot_breaks

def refine_perm_breaks_only(perm, gt, A_init, A_t, W_init, max_iters=2000):
    best_perm = perm.copy()
    best = eval_breaks_only(best_perm, A_init, A_t, W_init)
    best_obj = best[1]  # breaks_w
    iters = 0
    improved = True
    while improved and iters < max_iters:
        improved = False; iters += 1
        for z in sorted(gt.keys()):
            tgt_list = gt[z]
            L = len(tgt_list)
            if L < 2: continue
            for a in range(L):
                i = tgt_list[a]
                for b in range(a+1, L):
                    j = tgt_list[b]
                    best_perm[i], best_perm[j] = best_perm[j], best_perm[i]
                    k_w, br_w, k, br = eval_breaks_only(best_perm, A_init, A_t, W_init)
                    obj = br_w
                    if (obj + 1e-12 < best_obj) or (abs(obj - best_obj) < 1e-12 and k > best[2]):
                        best = (k_w, br_w, k, br); best_obj = obj
                        improved = True
                    else:
                        best_perm[i], best_perm[j] = best_perm[j], best_perm[i]
                if improved: break
            if improved: break
    return best_perm, *best

def refine_perm_breaks_across_images(perm, gt, A_init, W_init, A_list, max_iters=2000):
    best_perm = perm.copy()
    best = total_breaks_across_images(best_perm, A_init, W_init, A_list)
    best_obj = best[1]  # total breaks_w
    iters = 0
    improved = True
    while improved and iters < max_iters:
        improved = False; iters += 1
        for z in sorted(gt.keys()):
            tgt_list = gt[z]
            L = len(tgt_list)
            if L < 2: continue
            for a in range(L):
                i = tgt_list[a]
                for b in range(a+1, L):
                    j = tgt_list[b]
                    best_perm[i], best_perm[j] = best_perm[j], best_perm[i]
                    tot = total_breaks_across_images(best_perm, A_init, W_init, A_list)
                    obj = tot[1]
                    if obj + 1e-12 < best_obj or (abs(obj - best_obj) < 1e-12 and tot[2] > best[2]):
                        best = tot; best_obj = obj
                        improved = True
                    else:
                        best_perm[i], best_perm[j] = best_perm[j], best_perm[i]
                if improved: break
            if improved: break
    return best_perm, *best

# =======================
# Seed mapping (endpoints)
# =======================
def build_perm_seed_endpoints(mobile, target, cov_scale, base_mode, bond_weight_exp, pair_mult_sym, k_dist=4, local_refine=True, max_refine_iters=2000):
    Zm, Zt = mobile.numbers, target.numbers
    if sorted(Zm.tolist()) != sorted(Zt.tolist()):
        raise ValueError("Compositions differ; cannot build a permutation.")

    A_m, D_m = _adjacency_by_covalent(mobile, cov_scale)
    A_t, D_t = _adjacency_by_covalent(target, cov_scale)
    W_m = weight_matrix_depth2(mobile, A_m, base_mode, bond_weight_exp, pair_mult_sym)
    W_t = weight_matrix_depth2(target, A_t, base_mode, bond_weight_exp, pair_mult_sym)

    gm = defaultdict(list); gt = defaultdict(list)
    for i, z in enumerate(Zm): gm[int(z)].append(i)
    for j, z in enumerate(Zt): gt[int(z)].append(j)

    perm = np.empty(len(mobile), dtype=int)  # target idx -> mobile idx

    # Hungarian per element using endpoint costs
    for z in sorted(gt.keys()):
        idx_m = gm[z]; idx_t = gt[z]
        if len(idx_m) != len(idx_t):
            raise ValueError(f"Element Z={z} count mismatch.")
        M = len(idx_m)
        C = np.zeros((M, M), float)
        for a, i in enumerate(idx_m):
            for b, j in enumerate(idx_t):
                C[a, b] = _atom_cost_endpoint(i, j, A_m, A_t, D_m, D_t, W_m, W_t, Zm, Zt, k_dist=k_dist)
        r, c = linear_sum_assignment(C)
        for a, b in zip(r, c):
            perm[idx_t[b]] = idx_m[a]

    # Optional 2-swap refine vs final only
    keep_w, breaks_w, keep, br = eval_breaks_only(perm, A_m, A_t, W_m)
    if local_refine:
        perm, keep_w, breaks_w, keep, br = refine_perm_breaks_only(
            perm, gt, A_m, A_t, W_m, max_iters=max_refine_iters
        )
    return perm, keep_w, breaks_w, keep, br

# =======================
# NEB helpers
# =======================
def build_idpp_images(initial_aligned, final, n_images_half, parallel=False):
    n_images = 2 * n_images_half
    images = [initial_aligned] + [initial_aligned.copy() for _ in range(n_images)] + [final.copy()]
    neb = NEB(images, climb=True, parallel=parallel)
    neb.interpolate(method="idpp")
    return images, neb

def adjacencies_for_images(images, cov_scale):
    A_list, D_list = [], []
    for img in images:
        A, D = _adjacency_by_covalent(img, cov_scale)
        A_list.append(A); D_list.append(D)
    return A_list, D_list

# =======================
# Pipeline
# =======================
def run_neb_depth2_pipeline(
    initial_file,
    final_file,
    outdir,
    fairchem_model,
    device,
    cov_scale_init,
    cov_scale_inter,
    n_images_half,
    endpoint_relax_fmax,
    endpoint_relax_steps,
    neb_relax_fmax,
    neb_relax_steps,
    weight_base_mode,
    bond_weight_exp,
    depth1_coef,
    depth2_coef,
    pair_mult_sym,
    k_dist,
    max_mapping_iters,
    show_preview_before,
    show_preview_after,
    use_masses_kabsch,
    spin,
    charge,
    neb_parallel
):
    os.makedirs(outdir, exist_ok=True)

    initial = read(initial_file)
    final   = read(final_file)

    if spin is not None:  # optional metadata
        initial.info["spin"] = spin
        final.info["spin"] = spin
    if charge is not None:
        initial.info["charge"] = charge
        final.info["charge"]   = charge

    # Build initial graphs/weights (reference for breaks is INITIAL)
    A_init, _ = _adjacency_by_covalent(initial, cov_scale_init)
    # Inject custom depth coefficients via partial application style
    def _W(atoms, A):
        # temporarily patch global-like depth coefs by wrapping _depth2_node_score via args
        # simplest: pass depth1/2 via lambda closure in weight computation
        Z = atoms.numbers.astype(int)
        n = len(Z); W = np.zeros((n, n), float)
        for i in range(n):
            for j in np.where(A[i])[0]:
                if j <= i: continue
                si = _depth2_node_score(i, A, Z, weight_base_mode, exclude=j, d1c=depth1_coef, d2c=depth2_coef)
                sj = _depth2_node_score(j, A, Z, weight_base_mode, exclude=i, d1c=depth1_coef, d2c=depth2_coef)
                w = (si * sj) ** 0.5
                if bond_weight_exp != 1.0: w = w ** bond_weight_exp
                mult = pair_mult_sym.get(_pair_sym(Z[i], Z[j]), 1.0)
                W[i, j] = W[j, i] = w * mult
        return W

    W_init = _W(initial, A_init)

    # Seed mapping (endpoints)
    perm, keep_w, breaks_w, keep, br = build_perm_seed_endpoints(
        initial, final,
        cov_scale=cov_scale_init,
        base_mode=weight_base_mode,
        bond_weight_exp=bond_weight_exp,
        pair_mult_sym=pair_mult_sym,
        k_dist=k_dist,
        local_refine=True,
        max_refine_iters=2000
    )
    print(f"[seed] vs final: keep_w={keep_w:.3f} breaks_w={breaks_w:.3f} | keep={keep} breaks={br}")

    # Relax endpoints with FAIR-Chem
    predictor = make_fairchem_predictor(fairchem_model, device)
    for atoms, tag in [(initial, "initial"), (final, "final")]:
        attach_fairchem_calc(atoms, predictor)
        opt = LBFGS(atoms, logfile=os.path.join(outdir, f"{tag}_opt.log"))
        opt.run(fmax=endpoint_relax_fmax, steps=endpoint_relax_steps)
        write(os.path.join(outdir, f"{tag}_relaxed.xyz"), atoms)

    # Recompute initial graph after relaxation
    A_init, _ = _adjacency_by_covalent(initial, cov_scale_init)
    W_init = _W(initial, A_init)

    # Outer iteration: IDPP images + refine mapping across ALL images
    images = None
    neb = None

    for it in range(max_mapping_iters):
        print(f"\n=== Iteration {it} ===")

        initial_reindexed = initial[perm]
        write(os.path.join(outdir, f"initial_reindexed_iter{it}.xyz"), initial_reindexed)
        initial_aligned = kabsch_align(initial_reindexed, final, use_masses=use_masses_kabsch)
        write(os.path.join(outdir, f"initial_aligned_iter{it}.xyz"), initial_aligned)
        print(f"[iter {it}] RMSD = {rmsd(initial_aligned, final):.4f} Å")

        images, neb = build_idpp_images(initial_aligned, final, n_images_half, parallel=neb_parallel)

        if show_preview_before:
            try:
                # Avoid interactive viewer in headless; dump traj instead
                from ase.visualize import view
                view(images)
            except Exception:
                Trajectory(os.path.join(outdir, f"neb_preview_iter{it}.traj"), 'w').write_images(images)

        A_list, _ = adjacencies_for_images(images, cov_scale_inter)

        gt = defaultdict(list)
        for j, z in enumerate(images[0].numbers):
            gt[int(z)].append(j)

        gperm, tot_keep_w, tot_breaks_w, tot_keep, tot_breaks = refine_perm_breaks_across_images(
            perm, gt, A_init, W_init, A_list, max_iters=2000
        )

        map_path = os.path.join(outdir, f"mapping_iter{it}_target_to_initial.txt")
        with open(map_path, "w") as fh:
            fh.write("# target_index  initial_original_index\n")
            for t_idx, m_idx in enumerate(gperm.tolist()):
                fh.write(f"{t_idx:6d} {m_idx:6d}\n")

        print(f"[iter {it}] TOTAL keep_w={tot_keep_w:.3f} breaks_w={tot_breaks_w:.3f} | keep={tot_keep} breaks={tot_breaks}")

        if np.array_equal(gperm, perm):
            print(f"[iter {it}] permutation converged.")
            break
        perm = gperm

    # Attach FAIR-Chem to images and run NEB optimization
    for idx, img in enumerate(images):
        attach_fairchem_calc(img, predictor)
        write(os.path.join(outdir, f'neb_image_{idx:02d}_preopt.xyz'), img)

    neb_opt = LBFGS(neb, trajectory=os.path.join(outdir, "neb.traj"), logfile=os.path.join(outdir, "neb_opt.log"))
    neb_opt.run(fmax=neb_relax_fmax, steps=neb_relax_steps)

    energies = [img.get_potential_energy() for img in images]
    for idx, img in enumerate(images):
        write(os.path.join(outdir, f'neb_image_{idx:02d}.xyz'), img)

    if show_preview_after:
        try:
            from ase.visualize import view
            view(images)
        except Exception:
            Trajectory(os.path.join(outdir, "neb_final_preview.traj"), 'w').write_images(images)

    # Plot energy profile
    plt.figure(figsize=(8, 5))
    plt.plot(range(len(energies)), energies, 'o-', linewidth=2, markersize=8)
    plt.xlabel('Image Index', fontsize=14)
    plt.ylabel('Potential Energy (eV)', fontsize=14)
    plt.title('NEB Path Energy Profile', fontsize=16)
    plt.grid(True, linestyle='--', alpha=0.3)
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    plt.tight_layout()
    plot_path = os.path.join(outdir, "neb_energy_profile.png")
    plt.savefig(plot_path, dpi=200)

    i_ts = int(np.argmax(energies))
    fwd = energies[i_ts] - energies[0]
    rev = energies[i_ts] - energies[-1]
    print(f"[NEB] TS image ~{i_ts}, forward barrier={fwd:.4f} eV, reverse barrier={rev:.4f} eV")
    print(f"[NEB] Energies written. Plot: {plot_path}")

# =======================
# CLI
# =======================
def main():
    p = argparse.ArgumentParser(description="Depth-2-weighted NEB path with minimal bond breaks and FAIR-Chem relaxation")
    p.add_argument("initial_xyz")
    p.add_argument("final_xyz")
    p.add_argument("--outdir", default="neb_depth2", help="Output directory")
    p.add_argument("--fairchem-model", default="uma-s-1p1")
    p.add_argument("--device", choices=["cuda", "cpu"], default=("cuda" if torch.cuda.is_available() else "cpu"))

    # Endpoint + NEB settings
    p.add_argument("--cov-scale-init", type=float, default=1.20)
    p.add_argument("--cov-scale-inter", type=float, default=1.20)
    p.add_argument("--n-images-half", type=int, default=5)
    p.add_argument("--endpoint-relax-fmax", type=float, default=0.05)
    p.add_argument("--endpoint-relax-steps", type=int, default=100)
    p.add_argument("--neb-relax-fmax", type=float, default=0.05)
    p.add_argument("--neb-relax-steps", type=int, default=100)
    p.add_argument("--neb-parallel", action="store_true")

    # Depth-2 weights
    p.add_argument("--weight-base-mode", choices=["mass", "Z"], default="mass")
    p.add_argument("--bond-weight-exp", type=float, default=1.0)
    p.add_argument("--depth1-coef", type=float, default=1.0)
    p.add_argument("--depth2-coef", type=float, default=0.5)
    p.add_argument("--k-dist", type=int, default=4)

    # Outer loop
    p.add_argument("--max-mapping-iters", type=int, default=10)

    # Visualization toggles
    p.add_argument("--show-preview-before", action="store_true")
    p.add_argument("--show-preview-after", action="store_true")
    p.add_argument("--use-masses-kabsch", action="store_true", help="Mass-weighted alignment (default off)")

    # Optional metadata
    p.add_argument("--spin", type=int, default=None)
    p.add_argument("--charge", type=int, default=None)

    args = p.parse_args()

    # Pair multipliers placeholder (edit here if needed)
    pair_mult_sym = {
        # ('C','O'): 1.2,
        # ('C','H'): 0.8,
    }

    run_neb_depth2_pipeline(
        initial_file=args.initial_xyz,
        final_file=args.final_xyz,
        outdir=args.outdir,
        fairchem_model=args.fairchem-model if hasattr(args, "fairchem-model") else args.fairchem_model,
        device=args.device,
        cov_scale_init=args.cov_scale_init,
        cov_scale_inter=args.cov_scale_inter,
        n_images_half=args.n_images_half,
        endpoint_relax_fmax=args.endpoint_relax_fmax,
        endpoint_relax_steps=args.endpoint_relax_steps,
        neb_relax_fmax=args.neb_relax_fmax,
        neb_relax_steps=args.neb_relax_steps,
        weight_base_mode=args.weight_base_mode,
        bond_weight_exp=args.bond_weight_exp,
        depth1_coef=args.depth1_coef,
        depth2_coef=args.depth2_coef,
        pair_mult_sym=pair_mult_sym,
        k_dist=args.k_dist,
        max_mapping_iters=args.max_mapping_iters,
        show_preview_before=args.show_preview_before,
        show_preview_after=args.show_preview_after,
        use_masses_kabsch=args.use_masses_kabsch,
        spin=args.spin,
        charge=args.charge,
        neb_parallel=args.neb_parallel
    )

if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    main()
