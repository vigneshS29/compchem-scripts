import os, sys, argparse, subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()

parser.add_argument('-f',dest='file_path', type=str, default='data.csv',help='Input File, Default = data.csv')
parser.add_argument('-o',dest='output', type=str, default='energy_diagram.pdf',help='Output Image Name, Default = energy_diagram.pdf')
parser.add_argument('-t',dest='func', type=str, default='gaussian',help='Spline Function (gaussian/exponential), Default = gaussian')
args = parser.parse_args()

def main(argv):

    color_list = ['red','blue']

    data = read_data(args.file_path)

    plt.figure(figsize=(20,5))

    for j in data:
        x,y = data[j][0],data[j][1]

        for i in range(0,len(x)-1,1):
            xtemp = [x[i],x[i+1]]
            ytemp = [y[i],y[i+1]]
            
            if args.func == 'exponential': new_x,new_y = exponential_func(xtemp,ytemp)
            if args.func == 'gaussian': new_x,new_y = gaussian_func(xtemp,ytemp)
            
            plt.plot(new_x, new_y,color='black',alpha=0.5)
            plt.scatter(xtemp,ytemp,color=color_list[int(j)-1],edgecolor='black',s=50)
        
    plt.xticks(np.arange(0,max(x)+1))
    plt.xlim(min(x)-1,max(x)+1)
    plt.ylim(min(y)-5,max(y)+5)
    plt.xlabel('Reaction Coordinate',fontsize=15)
    plt.ylabel('Energy (kcal/mol)',fontsize=15)

    plt.savefig(args.output) 

    return 

def read_data(file_path):
    
    data = {}
    with open(file_path) as f:
        f.readline()
        for fields in f:
            i = fields.rstrip().split(',')
            
            if i[3] in data: 
                data[i[3]][0] += [float(i[0])]
                data[i[3]][1] += [float(i[1])]
            else:
                data[i[3]] = [[],[]]
                data[i[3]][0] += [float(i[0])]
                data[i[3]][1] += [float(i[1])]
    
    return data

def exponential_func(x,y):
    
    if y[1] > y[0]:
        xtemp = np.arange(min(x),max(x)+0.01,0.1)
        ytemp = np.exp(0.2*xtemp**2)
        ytemp = ytemp - min(ytemp)
        ytemp = ytemp/max(ytemp)
        ytemp = (max(y) - min(y))*ytemp + min(y)
    
    else:
        xtemp = np.arange(min(x),max(x)+0.01,0.1)
        ytemp = np.exp(-0.2*(xtemp**2))
        ytemp = ytemp - min(ytemp)
        ytemp = ytemp/max(ytemp)
        ytemp = (max(y) - min(y))*ytemp + min(y)
           
    return xtemp,ytemp

def gaussian_func(x,y):
    
    if y[1] > y[0]:
        sig = 2
        mu = max(x)
        xtemp = np.arange(min(x),max(x)+0.01,0.1)
        ytemp = np.exp(-((xtemp-mu*np.ones(len(xtemp)))/sig)**2)
        ytemp = ytemp - min(ytemp)
        ytemp = ytemp/max(ytemp)
        ytemp = (max(y) - min(y))*ytemp + min(y)
        
    else:
        sig = 2
        mu = max(x)
        xtemp = np.arange(min(x),max(x)+0.01,0.1)
        ytemp = np.exp(-((xtemp-mu*np.ones(len(xtemp)))/sig)**2)
        ytemp = ytemp - min(ytemp)
        ytemp = ytemp/max(ytemp)
        ytemp = (max(y) - min(y))*ytemp + min(y)
        y = ytemp[::-1]
        ytemp = y
        
    return xtemp,ytemp

if __name__ == "__main__":
    main(sys.argv[1:])
