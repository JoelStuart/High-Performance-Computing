#!/bin/bash -l
#SBATCH --job-name=SerProj
#SBATCH --partition=iceq
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20000
#SBATCH --time=00:43:00
#SBATCH --mail-type=FAIL
# ...end of SLURM preamble...
# ...load required modules...
module load intel intel-mpi
# ...now run your tasks. Don't forget srun!
echo "This is job '$SLURM_JOB_NAME' (id: $SLURM_JOB_ID) running on the following nodes:"
echo $SLURM_NODELIST
echo "Now we start the show:"
gcc -c projSer.c
gcc -o projSer projSer.o -I. -lm -O2
export TIMEFORMAT="%E sec"
time ./projSer

echo "All done."

