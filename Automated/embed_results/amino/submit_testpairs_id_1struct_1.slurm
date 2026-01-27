#!/bin/sh
#SBATCH --job-name=testpairs_id_1struct_1.inp
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=8
#SBATCH --mem=16000
#SBATCH --error=/lustre/hpc/kemi/hteahan/Orca_workflow/Automated/embed_results/amino/testpairs_id_1struct_1.inp.err
#SBATCH --output=/lustre/hpc/kemi/hteahan/Orca_workflow/Automated/embed_results/amino/testpairs_id_1struct_1.inp.out
#SBATCH --time=08:00:00
#SBATCH --partition=kemi1
#SBATCH --no-requeue
module load mpi/openmpi-x86_64

#Move to scratch dir for calculations
cp testpairs_id_1struct_1.inp /scratch/$SLURM_JOB_ID/
cd /scratch/$SLURM_JOB_ID/
#Run calculations
/groups/kemi/hteahan/opt/orca_6_1_0_linux_x86-64_shared_openmpi418/orca testpairs_id_1struct_1.inp

mv * /lustre/hpc/kemi/hteahan/Orca_workflow/Automated/embed_results/amino
