#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --time=00:35:59
#SBATCH --job-name=POST
#SBATCH --output=LOG/POST
#SBATCH --account=rrg-liuqy
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
