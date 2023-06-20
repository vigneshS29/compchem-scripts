#!/bin/python

#Author: Vignesh Sathyaseelan (vsathyas@purdue.edu)

import os, sys, argparse, subprocess
import multiprocessing as mp
from multiprocessing import Pool
from rdkit import Chem,RDLogger  
from rdkit.Chem import AllChem

parser = argparse.ArgumentParser()
parser.add_argument('-s',dest='smile', type=str, help= 'space delimited smile strings')
args = parser.parse_args()

def main(argv):

    arr = [(i,count_i) for count_i,i in enumerate(args.smile.split())]
    
    with Pool(mp.cpu_count()) as p:
        result = p.starmap(get_xyz_frm_smi,arr)

    return

def check_validity(smile):
    
    from rdkit import RDLogger  
    RDLogger.DisableLog('rdApp.*')
    validity = False
    m = Chem.MolFromSmiles(smile)
    if m != None:
        validity =  True
    return validity

def get_xyz_frm_smi(smile,out_name,out_dir='xyz_folder',steps=10000,forcefield='uff'):
    
    xyz_folder = "{}/{}".format(os.getcwd(),out_dir)
    
    if os.path.exists(xyz_folder) is False: 
        os.mkdir(xyz_folder)
    
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

    else: print(f'{smile} is a invalid smile')
        
    return 
    
if __name__ == '__main__':
    main(sys.argv[1:])
