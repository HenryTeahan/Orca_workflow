#!/bin/sh
#SBATCH --job-name=testpairs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --mem=20000
#SBATCH --error=/lustre/hpc/kemi/hteahan/Orca_workflow/Automated/embed_results/amino/testpairs.err
#SBATCH --output=/lustre/hpc/kemi/hteahan/Orca_workflow/Automated/embed_results/amino/testpairs.out
#SBATCH --time=08:00:00
#SBATCH --partition=kemi1
#SBATCH --no-requeue
#micromamba activate Complex
#Move to scratch dir for calculations
cp /groups/kemi/hteahan/Orca_workflow/Automated/embed_results/testpairs.csv /scratch/$SLURM_JOB_ID/
ln -s /lustre/hpc/kemi/hteahan/Orca_workflow/Automated/scripts/TMC_embed /scratch/$SLURM_JOB_ID/TMC_embed
cp /lustre/hpc/kemi/hteahan/Orca_workflow/Automated/scripts/embed.py /scratch/$SLURM_JOB_ID/
ln -s /lustre/hpc/kemi/hteahan/opt/xtb-6.7.1/xtb-dist/bin/xtb /scratch/$SLURM_JOB_ID/xtb
cd /scratch/$SLURM_JOB_ID/
export PYTHONPATH=/scratch/$SLURM_JOB_ID:$PYTHONPATH


#Run calculations
micromamba run -n Complex python embed.py --xyz_files None --smiles /groups/kemi/hteahan/Orca_workflow/Automated/embed_results/testpairs.csv --smiles_col smiles_0 --N_tries 20 --xtb_path /scratch/$SLURM_JOB_ID/xtb
cp * /lustre/hpc/kemi/hteahan/Orca_workflow/Automated/embed_results/amino
cp * /lustre/hpc/kemi/hteahan/Orca_workflow/Automated/embed_results/amino
