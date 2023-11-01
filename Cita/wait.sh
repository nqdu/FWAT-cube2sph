#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:16:00
#SBATCH --job-name=wait
#SBATCH --output=wait_%j.txt
#SBATCH --partition=compute


echo "waiting for successed"
