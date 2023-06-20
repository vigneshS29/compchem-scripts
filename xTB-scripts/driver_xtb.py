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

if __name__ == '__main__':
    main(sys.argv[1:])