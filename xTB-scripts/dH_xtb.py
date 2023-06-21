import os, sys, argparse, subprocess, pickle, numpy as np
from multiprocessing import Pool
from xtb_functions import *

parser = argparse.ArgumentParser()
parser.add_argument('-f',dest='in_file_name' ,type=str, default = 'test.p', help ='file to read')
args = parser.parse_args()


def main(argv):

    with open(args.in_file_name, 'rb') as f:
        
        try:
            os.mkdir('{}/{}'.format(os.getcwd(),args.in_file_name.split('/')[-1].split('.')[0]))
        except:
            pass

        os.chdir('{}/{}'.format(os.getcwd(),args.in_file_name.split('/')[-1].split('.')[0]))

        
        x = pickle.load(f)
        xyz_write(name='reactant.xyz',elements=np.array(x[list(x.keys())[0]]['prop']['E']),geo=np.array([x[list(x.keys())[0]]['prop']['G']]))
        E0,_,_ = xtb_geo_opt('reactant.xyz')
        os.system('rm *.xyz')
        E = []

        count = 0
        for count_i,i in enumerate(list(x[list(x.keys())[0]]['possible_products'].keys())):

            xyz_write(name='{}.xyz'.format(i),elements=np.array(x[list(x.keys())[0]]['possible_products'][i]['E']),geo=np.array(x[list(x.keys())[0]]['possible_products'][i]['G_list']))
            E1,_,_ = xtb_geo_opt('{}.xyz'.format(i))
            if (E1-E0)*627.5 < 0.0: count +=1 
            E+= [float((E1-E0)*627.5)]
            os.system('rm *.xyz')

        print(args.in_file_name.split('/')[-1].split('.')[0],count,np.sort(np.array(E))[:5])


    return

def xyz_write(name,elements,geo,append_opt=False,comment=''):

    if append_opt == True:
        open_cond = 'a'
    else:
        open_cond = 'w'
    
    geo = geo[0]
    
    with open(name,open_cond) as f:
        f.write('{}\n'.format(len(elements)))
        f.write('{}\n'.format(comment))
        for count_i,i in enumerate(elements):
            f.write("{:<20s} {:< 20.8f} {:< 20.8f} {:< 20.8f}\n".format(i,geo[count_i][0],geo[count_i][1],geo[count_i][2]))
    return 

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
