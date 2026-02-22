#!/bin/bash

#SBATCH -J at12
#SBATCH -n 1 # Number of cores requested
#SBATCH -t 2:00:00 # Runtime
#SBATCH -p conroy_priority # Partition to submit to
#SBATCH --constraint="intel"
#SBATCH --mem=4000 # Memory MB
#SBATCH -o run.o # Standard out goes to this file
#SBATCH -e run.e # Standard err goes to this file

cd $SCRATCH/conroy_lab/cconroy/kurucz/atlas12

#atlas12.sh start end element abundance_change output_dir input_dir
./atlas12.sh $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_ID 1 0.00 c3k_v2.2/at12_feh+0.00_afe+0.0 c3k_v2.1/at12_feh+0.00_afe+0.0 >& out.$SLURM_ARRAY_TASK_ID

