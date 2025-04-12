#!/bin/bash

job2=$(sbatch tmp.fwat2.sh | cut -d ' ' -f4)
job3=$(sbatch --dependency=afterok:${job2} tmp.fwat3.tele.sh | cut -d ' ' -f4)
job4=$(sbatch --dependency=afterok:${job2} tmp.fwat3.noise.sh | cut -d ' ' -f4)
sbatch --dependency=afterok:${job3},${job4} tmp.fwat4.sh
