

#!/bin/bash
iter=1
iterMax=20
myExe=ctmo
nprocess=2

rm tktilde.arma tloc.arma hybFM.arma config*.dat
while [  $iter -lt $iterMax ];
    do
        mpirun -np $nprocess $myExe params $iter
	    iter=$[iter+1]
    done


