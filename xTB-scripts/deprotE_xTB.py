import numpy as np, copy, os, sys, subprocess, shutil, time
from rdkit import AllChem, Chem, DataStructs

def main(argv):

    smile_l = ['FC(F)(F)COCC(F)(F)F', 'COCCC#N', 'C1=CC(=CC(=C1)F)C#N','C1=CC(=NC(=C1)F)F','COCC#N', 'CCC#N', 'C1=CC=C(C(=C1)C#N)F', 
           'CC(=O)N(C)C', 'C1=CC(=C(C(=C1)F)F)C#N','C','CC#N']

    for smile in smile_l:
        print(smile,calc_deprotE(smile,deprotonate(smile)))

    return

def calc_deprotE(ref_mol,deprot_mols):
    

    smile_file = ref_mol.replace('(','_').replace(')','_')
    get_xyz_frm_smi('{}'.format(ref_mol),out_dir='xyz_folder_deprotE',out_name=f'{smile_file}',charge=0)
    E0,_,_ = xtb_geo_opt(f'xyz_folder_deprotE/{smile_file}.xyz',charge=0)

    smile_file = '[H]'.replace('(','_').replace(')','_')
    get_xyz_frm_smi('{}'.format('[H]'),out_dir='xyz_folder_deprotE',out_name=f'{smile_file}',charge=0)
    E1,_,_ = xtb_geo_opt(f'xyz_folder_deprotE/{smile_file}.xyz',charge=1)

    E = []
    for count_i,i in enumerate(deprot_mols):
        smile_file = i.replace('(','_').replace(')','_')
        get_xyz_frm_smi('{}'.format(i),out_dir='xyz_folder_deprotE',out_name=f'{smile_file}',charge=0) 
        E2,_,_ = xtb_geo_opt(f'xyz_folder_deprotE/{smile_file}.xyz',charge=-1)
    
        E+= [float((E0 - E2 - E1))]
    

    return np.sort(np.array(E))[0]

def deprotonate(smile):
    mol_list = []
    mol_ref = copy.deepcopy(Chem.MolFromSmiles(smile))
    mol = copy.deepcopy(Chem.MolFromSmiles(smile))
    for atom_count, atom in enumerate(mol.GetAtoms()):
            if mol.GetAtomWithIdx(atom_count).GetTotalNumHs() > 0:
                    mol.GetAtomWithIdx(atom_count).SetFormalCharge(-1)
            mol_list += [Chem.MolToSmiles(mol)]
            mol = Chem.MolFromSmiles(smile)
    mol_list = [_ for _ in mol_list if _ != Chem.MolToSmiles(Chem.MolFromSmiles(smile))]
   
    return list(set(mol_list))

def check_validity(smile):
    
    from rdkit import RDLogger  
    RDLogger.DisableLog('rdApp.*')
    validity = False
    m = Chem.MolFromSmiles(smile)
    if m != None:
        validity =  True
    return validity

def get_xyz_frm_smi(smile,out_dir='xyz_folder',out_name='0',steps = 10000,forcefield = 'uff',charge=0):
    xyz_folder = "{}/{}".format(os.getcwd(),out_dir)
    
    try:
        os.mkdir(xyz_folder)
    except:
        pass

    m = Chem.MolFromSmiles(smile)
    if check_validity(smile):
        m2= Chem.AddHs(m)
        AllChem.EmbedMolecule(m2)
        tmp_filename = f'{out_name}.tmp.mol'
        with open(tmp_filename,'w') as g: g.write(Chem.MolToMolBlock(m2))
        xyzf = '{}/{}.xyz'.format(xyz_folder,out_name)
        cmd = 'obabel {} -O {} --sd --minimize --steps {} --ff {} --gen3d'.format(tmp_filename,xyzf,steps,forcefield)
        output = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,bufsize=-1).wait()
        os.remove(tmp_filename)

        if charge != 0:
            with open(xyzf,'r') as f:
                data = f.readlines()
            
            data[1] = 'q {}\n'.format(charge)
        
            os.remove(xyzf)

            with open(xyzf,'w') as f:
                for i in data:
                    f.write(i)
        
        return True


def xtb_geo_opt(xyz_file,charge=-1,unpair=0,niter=100,accuracy=1.0,namespace='xTB',workdir='.',level='normal',fixed_atoms=[],output_xyz=None,cleanup=False):

    current_dir = os.getcwd()

    # change to xTB working directory
    os.chdir(workdir)

    # write input file is needed
    if len(fixed_atoms) > 0:
        with open('{}-xtb.inp'.format(namespace),'w') as f:
            f.write('$fix\n')
            for ind in fixed_atoms:
                f.write('atoms: {}\n'.format(ind))
            f.write('$end\n')

        # generate command line
        substring= "xtb -c {} -u {} -a {} --input {}-xtb.inp --iterations {} --opt {} --namespace '{}' {}"
        code_exe = substring.format(charge,unpair,accuracy,namespace,niter,level,namespace,xyz_file)
        
    else:
        # generate command line
        substring= "xtb -c {} -u {} -a {} --iterations {} --opt {} --namespace '{}' {}" 
        code_exe = substring.format(charge,unpair,accuracy,niter,level,namespace,xyz_file)

    # run xtb
    output = subprocess.Popen(code_exe, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]
    output = output.decode('utf-8')
    Energy = 0.0
    with open("{}_xTB-opt.txt".format(namespace),"w") as g: g.write(output)
        
    with open("{}_xTB-opt.txt".format(namespace),"r") as g:
        for lc,lines in enumerate(g):
            fields=lines.split()
            if len(fields)==6 and fields[1] =="TOTAL" and fields[2] == "ENERGY":
                Energy = float(fields[3])

    opt_xyz_file = "{}/{}.xtbopt.xyz".format(workdir,namespace)
    if Energy != 0.0:
    
        # clean up the xtb files is the flag is on
        if cleanup:
            if os.path.isfile("{}/{}.xtbopt.xyz".format(current_dir,namespace)) is False:
                shutil.move(opt_xyz_file,current_dir)
            subprocess.Popen("rm {}.*".format(namespace),shell=True)            
            time.sleep(0.1)
            if os.path.isfile("{}/{}.xtbopt.xyz".format(workdir,namespace)) is False:
                shutil.move("{}/{}.xtbopt.xyz".format(current_dir,namespace),workdir)

        # change back to original folder
        os.chdir(current_dir)

        # copy the optput xyz file to given path if is needed
        if output_xyz is not None:
            shutil.copy2(opt_xyz_file,output_xyz)
            #print("Geometry optimization is done at xtb level with single point energy:{} and resulting xyz file {}".format(Energy,output_xyz))
            return Energy,output_xyz,True

        else:
            #print("Geometry optimization is done at xtb level with single point energy:{} and resulting xyz file {}".format(Energy,opt_xyz_file))
            return Energy,opt_xyz_file,True
    else:
        #print("xTB Geo-opt fails")
        # change back to original folder
        os.chdir(current_dir)
        return Energy,opt_xyz_file,False

if __name__ == '__main__':
    main(sys.argv[1:])