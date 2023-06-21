#!/bin/bash                                                                                 

#SBATCH --job-name $1
#SBATCH -N 1
#SBATCH -n $2
#SBATCH -t $3
#SBATCH -A $4
#SBATCH -o output.txt
#SBATCH -e error.txt

for file in $1*; do
    #will run several instances parallely 
    python /home/vsathyas/bin/GA/GA/driver_xtb.py -f /home/vsathyas/bin/GA/GA/$file >> data.txt & 
done

wait 
