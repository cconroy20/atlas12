#!/bin/bash

#SBATCH -J at12
#SBATCH -n 1 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 04:00:00 # Runtime
#SBATCH -p conroy_priority
#SBATCH --mem-per-cpu=4000 # Memory per cpu in MB (see also --mem)
#SBATCH -o /dev/null # Standard out goes to this file
#SBATCH -e /dev/null # Standard err goes to this file

cd workdir

../scripts/atlas12.sh 5 5 1 0.00 new test >& out




