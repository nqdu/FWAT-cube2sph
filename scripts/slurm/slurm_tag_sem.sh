#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=40
#SBATCH --array=1-4%5
#SBATCH --time=00:35:59
#SBATCH --job-name=FWD
#SBATCH --output=LOG/FWD-%j_set%a.txt
#SBATCH --account=rrg-liuqy
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=nanqiao.du@mail.utoronto.ca
