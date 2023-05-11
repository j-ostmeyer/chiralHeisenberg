#!/bin/bash
#SBATCH --job-name="heis_40"
#SBATCH --time=36:00:00
#SBATCH --exclusive=user
#SBATCH --account=DIRAC-DP137-CPU
#SBATCH --partition=skylake
#SBATCH -n 1 -c 4
#SBATCH --dependency=singleton

# Set up the environment
. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load slurm rhel7/global rhel7/default-peta4 intel/bundles/complib/2018.4

export OMP_NUM_THREADS=4
export OMP_PROC_BIND=close
export OMP_PLACES=cores

time ./xcheis >> $SLURM_JOB_ID.out
echo "numbers of processors: $OMP_NUM_THREADS" >> output
echo "Slurm Job ID: $SLURM_JOB_ID" >> output

echo $(date +"%s") end $SLURM_JOB_ID $(cat fort.100  2>/dev/null| wc -l ) \

