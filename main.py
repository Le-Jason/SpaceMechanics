from orbitInitPropagate import orbitInitPropagate
from orbitInitPropagate import perturbations
from Extra import multiPlot
from Extra import giveOrbitalElements
from Extra import readData
import numpy as np

def main():
    [Name,mass,radius,rotation,distance,eccentricity,J2] = readData('Earth')
    G = 0.0000000000000000000667430
    mu = G*mass
    state = [414+6371,0.0006189,51.6393,105.6372,234.1955,0]
    perts=perturbations()
    perts['J2'] + True
    Orbit1 = orbitInitPropagate('Earth',kepler=True,stateVec=state,perturbation=perts)
    t0 = 0
    y0 = [0,0,0,0,0,0]
    dt = 100
    tf = 100*5000
    tSol,ySol = Orbit1.propagateOrbit(t0,y0,dt,tf)
    Orbit1.plotKepler()


if __name__ == '__main__':
    main()