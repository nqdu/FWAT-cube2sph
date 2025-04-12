#!/bin/bash

job2=$(sbatch sbash_fwat2_postproc_opt.sh | cut -d ' ' -f4)
job3=$(sbatch --dependency=afterok:${job2} tmp.fwat3.tele.sh | cut -d ' ' -f4)
sbatch --dependency=afterok:${job3} sbash_fwat4_opt_model.sh
