#!/bin/bash
#SBATCH --job-name=run_workflow
#SBATCH --output=run_workflow.out
#SBATCH --error=run_workflow.err
#SBATCH --time=3-00:00:00        
#SBATCH --cpus-per-task=8        
#SBATCH --mem=8000               
#SBATCH --partition=kemi1        
python /groups/kemi/hteahan/Orca_workflow/Automated/scripts/run_job.py --DB_PATH /groups/kemi/hteahan/Orca_workflow/Automated/isomer0/db/jobs.db --SMILES /groups/kemi/hteahan/Orca_workflow/Automated/csvs/test.csv --SMILES_COL smiles_0 --CHARGE_COL charge_0 --ID_COL ligand_ID --JOB_ID test_job --SCRATCH_DIR /scratch/$SLURM_JOB_ID --xTB_path /groups/kemi/hteahan/opt/xtb-6.7.1/xtb-dist/bin/xtb

