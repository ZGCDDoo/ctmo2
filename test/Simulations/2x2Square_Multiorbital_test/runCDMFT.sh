

#!/bin/bash
iter=10
iterMax=20
myExe=../../ctmo

rm tktilde.arma tloc.arma hybFM.arma config*.dat
while [  $iter -lt $iterMax ];
    do
        $myExe params $iter
	    iter=$[iter+1]
    done


