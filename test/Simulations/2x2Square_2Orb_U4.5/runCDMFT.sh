

#!/bin/bash
ITER=1
ITERMAX=3
myExe=ctmo
nprocess=2

rm tktilde.arma tloc.arma hybFM.arma config*.dat
while [  $ITER -lt $ITERMAX ];
    do
        mpirun -np $nprocess $myExe params${ITER}.json
	    ITER=$[ITER+1]
    done


