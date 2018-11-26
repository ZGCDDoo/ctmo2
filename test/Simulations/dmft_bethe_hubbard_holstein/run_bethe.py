import subprocess as sp
import sys
import numpy as np
from shutil import copyfile


iter_start = int(sys.argv[1])
EXE = "mpirun -np 4 ctmo_nosc "
ww = 0.2

for ii in range(iter_start, iter_start + 30):
    cmd = EXE + " params " + str(ii)
    sp.run(cmd, shell=True)
    copyfile("params" + str(ii) + ".json", "params" + str(ii+1) + ".json")
    copyfile("greenUp.dat", "greenUp" + str(ii) + ".dat")
    copyfile("Obs.json", "Obs" + str(ii) + ".json")
    gf = np.loadtxt("greenUp.dat")
    gf_previous = np.loadtxt("greenUp" + str(ii-1) + ".dat")
    gf_weighted = (1.0 - ww)*gf + ww*gf_previous
    np.savetxt("greenUp.dat", gf_weighted)
