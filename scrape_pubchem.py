import argparse, sys
from pubchempy import Compound
from pubchemprops.pubchemprops import get_second_layer_props

parser = argparse.ArgumentParser()
parser.add_argument('-N',dest='num', type=int, default=1000, help= 'Numer of molecules to scrape')
parser.add_argument('-o',dest='outname', type=str, default='Tb', help= 'name of output file')
args = parser.parse_args()


def main(argv):

    for i in range(args.num):  
        
        with open('{}.csv'.format(args.outname),'w') as f:
            f.write('smile,Tb')
            try:  
                comp = Compound.from_cid(i)
                field = get_second_layer_props(comp.iupac_name,\
                                                ['Boiling Point'])['Boiling Point'][0]\
                                                    ['Value']['StringWithMarkup'][0]['String']

                if '°F' in field: T = (5/9)*(float(field.split()[0])-32)+273
                if '°C' in field: T = float(field.split()[0])+273

                f.write('{},{:<3.4}'.format(comp.isomeric_smiles,T))

            except: pass

    return

if __name__ == '__main__':
    main(sys.argv[1:])

