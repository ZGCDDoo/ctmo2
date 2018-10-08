

#!/bin/bash
iter=40
iterMax=50
myExe=../../ctmo

rm tktilde.arma tloc.arma hybFM.arma config*.dat
while [  $iter -lt $iterMax ];
    do
        $myExe params $iter
	    iter=$[iter+1]
    done


