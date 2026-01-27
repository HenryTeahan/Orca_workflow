#!/bin/bash
#SBATCH --job-name=run_workflow
#SBATCH --output=run_workflow.out
#SBATCH --error=run_workflow.err
#SBATCH --time=3-00:00:00        # 3 days max
#SBATCH --cpus-per-task=8        # lightweight, only DB and submission
#SBATCH --mem=8000                 # enough for Python
#SBATCH --partition=kemi1        # or your default partition

# Run the workflow manager
python /groups/kemi/hteahan/Orca_workflow/Automated/scripts/run_job.py --DB_PATH /groups/kemi/hteahan/Orca_workflow/Automated/db/jobs.db --SMILES /groups/kemi/hteahan/Orca_workflow/Automated/csvs/test.csv --SMILES_COL smiles_0 --CHARGE_COL charge_0 --ID_COL ligand_ID --JOB_ID test_job --SCRATCH_DIR /scratch/$SLURM_JOB_ID --xTB_path /groups/kemi/hteahan/opt/xtb-6.7.1/xtb-dist/bin/xtb

