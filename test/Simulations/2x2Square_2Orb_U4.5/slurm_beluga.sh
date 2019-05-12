#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --account=def-tremblay

module reset
module load nixpkgs/16.09 gcc/7.3.0 cmake armadillo/7.950.1 snappy boost-mpi cmake

ITER=1
ITERMAX=3
myExe=ctmo

if [ -a logfile ]
  then rm logfile
fi
rm tktilde*.arma tloc*.arma hybFM*.arma config*.dat

while [ $ITER -le $ITERMAX ]
do
  echo begin iteration $ITER at: `date` >> logfile 

  srun $myExe params${ITER}.json

  echo end iteration $ITER at: `date` >> logfile
  ITER=$[$ITER+1]
done

