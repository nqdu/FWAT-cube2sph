#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:16:00
#SBATCH --job-name=wait
#SBATCH --output=wait_%j.txt
#SBATCH --partition=compute


echo "waiting for success"
# for f in optimize/*step*;
# do
# 	rm -rf $f 
# done

# for f in solver/*step*;
# do
# 	rm -rf $f 
# done
