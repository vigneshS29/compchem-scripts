import sys,os,argparse,subprocess,shutil,time,glob,fnmatch,rdkit
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def main(argv):
    
    parser = argparse.ArgumentParser(description='This script takes a folder of reactant/product alignments and generate conformations .')
    parser.add_argument('coord_files', help = 'The program performs on a folder of given xyz file of reactants and products')
    parser.add_argument('-w', dest='working_folder', default='scratch',help = 'Working folder for performing xTB calculations')
    parser.add_argument('-N', dest='N_max', default=3, help = 'maximum number of conformation')
    args=parser.parse_args()

    N_max = int(args.N_max)

    input_folder = os.getcwd()+'/'+args.coord_files 

    if args.working_folder[0] != '/':
        work_folder = os.getcwd()+'/'+args.working_folder
    else:
        work_folder = args.working_folder

    if os.path.isdir(work_folder) is False: os.mkdir(work_folder)
    if os.path.isdir(work_folder+'/opt-folder') is False: os.mkdir(work_folder+'/opt-folder')
    if os.path.isdir(work_folder+'/ini-inputs') is False: os.mkdir(work_folder+'/ini-inputs')
    if os.path.isdir(work_folder+'/conf-folder') is False: os.mkdir(work_folder+'/conf-folder')

    input_files = sorted([os.path.join(dp, f) for dp, dn, filenames in os.walk(input_folder) for f in filenames if (fnmatch.fnmatch(f,"*.xyz") )])

    for count_k,k in enumerate(input_files):

        E,RG,PG     = parse_input(k)

        R_index     = k.split('/')[-1].split('.xyz')[0]+'-R'
        P_index     = k.split('/')[-1].split('.xyz')[0]+'-P'
        
        reactant_file = work_folder+'/opt-folder/'+'{}.xyz'.format(R_index)
        xyz_write(reactant_file, E, RG)

        product_file = work_folder+'/opt-folder/'+'{}.xyz'.format(P_index)
        xyz_write(product_file, E, PG)

        Rsmile = get_smi_frm_xyz(reactant_file)
        Psmile = get_smi_frm_xyz(product_file)

        R_m = Chem.AddHs(Chem.MolFromSmiles(Rsmile))
        P_m = Chem.AddHs(Chem.MolFromSmiles(Psmile))

        ps = AllChem.ETKDGv3()
        ps.randomSeed = 0xf00d

        R_cids = AllChem.EmbedMultipleConfs(R_m,N_max,params = ps)   
        P_cids = AllChem.EmbedMultipleConfs(P_m,N_max,params = ps)
        
        reactant_conf = [i for i in R_m.GetConformers()]
        product_conf = [i for i in P_m.GetConformers()]

        count = 0
        for count_Rconf,Rconf in enumerate(reactant_conf):

            RE = [R_m.GetAtomWithIdx(atom).GetSymbol() for atom in range(R_m.GetNumAtoms())]
            RG = np.array([[Rconf.GetAtomPosition(atom).x, Rconf.GetAtomPosition(atom).y, Rconf.GetAtomPosition(atom).z] for atom in range(R_m.GetNumAtoms())])
            xyz_write(work_folder+'/ini-inputs/'+'{}_{}.xyz'.format(R_index,count_Rconf),RE,RG)
            
            for count_Pconf,Pconf in enumerate(product_conf):
                
                PE = [P_m.GetAtomWithIdx(atom).GetSymbol() for atom in range(P_m.GetNumAtoms())]
                PG = np.array([[Pconf.GetAtomPosition(atom).x, Pconf.GetAtomPosition(atom).y, Pconf.GetAtomPosition(atom).z] for atom in range(P_m.GetNumAtoms())])
                xyz_write(work_folder+'/ini-inputs/'+'{}_{}.xyz'.format(P_index,count_Pconf),PE,PG)
                
                os.system('cat {} {} >> {}'.format(work_folder+'/ini-inputs/'+'{}_{}.xyz'.format(R_index,count_Rconf),work_folder+'/ini-inputs/'+'{}_{}.xyz'.format(P_index,count_Pconf),work_folder+'/conf-folder/'+'{}_{}.xyz'.format(k.split('/')[-1].split('.xyz')[0],count)))
                count+=1
        
    return

def get_smi_frm_xyz(xyz_file):
    cmd = 'obabel {} -osmi'.format(xyz_file)
    output = str(subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,bufsize=-1).communicate()[0].decode('utf8')).split()[0]
    return output

def xyz_write(name,elements,geo,append_opt=False,comment=''):

    if append_opt == True:
        open_cond = 'a'
    else:
        open_cond = 'w'
        
    with open(name,open_cond) as f:
        f.write('{}\n'.format(len(elements)))
        f.write('{}\n'.format(comment))
        for count_i,i in enumerate(elements):
            f.write("{:<20s} {:< 20.8f} {:< 20.8f} {:< 20.8f}\n".format(i,geo[count_i][0],geo[count_i][1],geo[count_i][2]))
    return 

def parse_input(input_xyz):

    name = input_xyz.split('/')[-1].split('xyz')[0]
    xyz  = ['','']
    count= 0

    # read in pairs of xyz file
    with open(input_xyz,"r") as f:
        for lc,lines in enumerate(f):
            fields = lines.split()
            if lc == 0: 
                N = int(fields[0])
                xyz[0] += lines
                continue

            if len(fields) == 1 and float(fields[0]) == float(N):
                count+=1

            xyz[count]+=lines

    with open('{}_reactant.xyz'.format(name),"w") as f:
        f.write(xyz[0])

    with open('{}_product.xyz'.format(name),"w") as f:
        f.write(xyz[1])

    # parse reactant info
    E,RG   = xyz_parse('{}_reactant.xyz'.format(name))

    # parse product info
    _,PG   = xyz_parse('{}_product.xyz'.format(name))
                
    try:
        os.remove('{}_reactant.xyz'.format(name))
        os.remove('{}_product.xyz'.format(name))
    except:
        pass

    return E,RG,PG

def xyz_parse(input):
    
    # Iterate over the remainder of contents and read the
    # geometry and elements into variable. Note that the
    # first two lines are considered a header
    with open(input,'r') as f:
        for lc,lines in enumerate(f):
            fields=lines.split()

            # Parse header
            if lc == 0:
                if len(fields) < 1:
                    print("ERROR in xyz_parse: {} is missing atom number information".format(input))
                    quit()
                else:
                    N_atoms = int(fields[0])
                    Elements = ["X"]*N_atoms
                    Geometry = np.zeros([N_atoms,3])
                    count = 0

            # Parse body
            if lc > 1:

                # Skip empty lines
                if len(fields) == 0:
                    continue            

                # Write geometry containing lines to variable
                if len(fields) > 3:

                    # Consistency check
                    if count == N_atoms:
                        print("ERROR in xyz_parse: {} has more coordinates than indicated by the header.".format(input))
                        quit()

                    # Parse commands
                    else:
                        Elements[count]=fields[0]
                        Geometry[count,:]=np.array([float(fields[1]),float(fields[2]),float(fields[3])])
                        count = count + 1

        # Consistency check
        if count != len(Elements):
            print("ERROR in xyz_parse: {} has less coordinates than indicated by the header.".format(input))

    return Elements,Geometry

if __name__ == "__main__":
    main(sys.argv[1:])
