from orbitInitPropagate import orbitInitPropagate
import numpy as np


r = 1000 + 6371
v = np.sqrt((0.0000000000000000000667430*5972000000000000000000000)/r)

r0 = [r,0,0]
v0 = [0,v,0]
y0 = r0 + v0
t0 = 0
N = 100*100
dt = 100


earthOrbit = orbitInitPropagate('Earth')
[tSol,ySol] = earthOrbit.propagateOrbit(t0,y0,dt,N)
earthOrbit.plot()



