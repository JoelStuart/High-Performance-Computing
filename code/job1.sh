#!/bin/bash -l
#SBATCH --job-name=ParaProj
#SBATCH --partition=iceq
#SBATCH --ntasks-per-node=4
#SBATCH --nodes=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=20000
#SBATCH --time=00:10:00
#SBATCH --mail-type=FAIL
# ...end of SLURM preamble...
# ...load required modules...
module load intel intel-mpi
# ...now run your tasks. Don't forget srun!
echo "This is job '$SLURM_JOB_NAME' (id: $SLURM_JOB_ID) running on the following nodes:"
echo $SLURM_NODELIST
echo "Now we start the show:"
make
export TIMEFORMAT="%E sec"
time mpirun -n 4 ./proj
time mpirun -n 8 ./proj


echo "All done."

