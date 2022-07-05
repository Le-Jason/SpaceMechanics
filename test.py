from orbitInitPropagate import orbitInitPropagate
from Extra import multiPlot
import numpy as np

def main():
    a = (409+418)/2 + 6371+1000
    e = 0.1
    i = 80
    RAAN = 45
    w = 45
    vu = 0

    stateVec = [a,e,i,RAAN,w,vu]
    ISSOrbit = orbitInitPropagate('Earth',kepler=True,stateVec=stateVec)
    t0 = 0
    y0 = [0,0,0,0,0,0]
    dt = 100
    tf = 100*500
    tSol,ySol = ISSOrbit.propagateOrbit(t0,y0,dt,tf)
    ISSOrbit.plot()

if __name__ == '__main__':
    main()