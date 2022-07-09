from orbitInitPropagate import orbitInitPropagate
from Extra import multiPlot
import numpy as np

def main():
    
    ISSOrbit = orbitInitPropagate('Earth',TLE=True,sat="GPS BIIR-2  (PRN 13)")
    t0 = 0
    y0 = [0,0,0,0,0,0]
    dt = 100
    tf = 100*500
    tSol,ySol = ISSOrbit.propagateOrbit(t0,y0,dt,tf)
    y = ySol
    ISSOrbit = orbitInitPropagate('Earth',TLE=True,sat="GPS BIIR-4  (PRN 20)")
    t0 = 0
    y0 = [0,0,0,0,0,0]
    dt = 100
    tf = 100*500
    tSol,ySol = ISSOrbit.propagateOrbit(t0,y0,dt,tf)
    ySol = [y,ySol] 
    Labels = ['GPS BIIR-2','GPS BIIR-4']
    Title = 'GPS'







    multiPlot('Earth',ySol,Labels,Title)
if __name__ == '__main__':
    main()