import numpy as np

exec(open("C:/Users/spong/Documents/Code/projects/spacemechanics/orbitInitPropagate.py").read())
exec(open("C:/Users/spong/Documents/Code/projects/spacemechanics/Extra.py").read())


perts=perturbations()
perts['J2'] + True
ISSOrbit = orbitInitPropagate('Earth',TLE=True,sat="GPS BIIR-2  (PRN 13)",perturbation=perts)
t0 = 0
y0 = [0,0,0,0,0,0]
dt = 100
tf = 100*500000
tSol,ySol = ISSOrbit.propagateOrbit(t0,y0,dt,tf)
y = ySol
ISSOrbit = orbitInitPropagate('Earth',TLE=True,sat="GPS BIIR-4  (PRN 20)")
tSol,ySol = ISSOrbit.propagateOrbit(t0,y0,dt,tf)
ySol = [y,ySol] 
Labels = ['GPS BIIR-2','GPS BIIR-4']
Title = 'GPS'

multiPlot('Earth',ySol,Labels,Title)