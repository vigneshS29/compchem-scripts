#!/bin/bash                                                                                 

#SBATCH --job-name $1
#SBATCH -n $2
#SBATCH -t 00:30:00
#SBATCH -A $3
#SBATCH -o output.txt
#SBATCH -e error.txt

for file in $1*; do
    #will run several instances parallely 
    python /home/vsathyas/bin/GA/GA/driver_xtb.py -f /home/vsathyas/bin/GA/GA/$file >> data.txt & 
done

wait 
