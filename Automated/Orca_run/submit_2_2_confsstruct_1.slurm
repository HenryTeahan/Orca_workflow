#!/bin/bash
#SBATCH --job-name=2_2_confsstruct_1.inp
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem 16000
#SBATCH --time=02:00:00
#SBATCH --output=/lustre/hpc/kemi/hteahan/Orca_workflow/Automated/Orca_run/2_2_confsstruct_1.inp.out
#SBATCH --error=/lustre/hpc/kemi/hteahan/Orca_workflow/Automated/Orca_run/2_2_confsstruct_1.inp.err
#SBATCH --partition=kemi1
#Move to scratch dir for calculations
module load mpi/openmpi-x86_64
cp 2_2_confsstruct_1.inp /scratch/$SLURM_JOB_ID/
cd /scratch/$SLURM_JOB_ID/

/groups/kemi/hteahan/opt/orca_6_1_0_linux_x86-64_shared_openmpi418/orca 2_2_confsstruct_1.inp

mv *.out /lustre/hpc/kemi/hteahan/Orca_workflow/Automated/Orca_run
mv *.err /lustre/hpc/kemi/hteahan/Orca_workflow/Automated/Orca_run
