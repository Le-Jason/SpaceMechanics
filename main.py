from orbitInitPropagate import orbitInitPropagate
from Extra import multiPlot
import numpy as np

def main():

    r = 10000 + 6371
    v = np.sqrt((0.0000000000000000000667430*5972000000000000000000000)/r)

    r0 = [r,0,0]
    v0 = [0,v,0]
    y0 = r0 + v0
    t0 = 0
    N = 100*500
    dt = 100


    earthOrbit = orbitInitPropagate('Earth')
    [tSol,ySol] = earthOrbit.propagateOrbit(t0,y0,dt,N)
    # earthOrbit.plot()
    r = 1000 + 6371
    v = np.sqrt((0.0000000000000000000667430*5972000000000000000000000)/r)+2
    r0 = [r,0,0]
    v0 = [0,v,0]
    print(v)
    y0 = r0 + v0
    earthOrbit2 = orbitInitPropagate('Earth')
    [tSol2,ySol2] = earthOrbit2.propagateOrbit(t0,y0,dt,N)
    ys = [ySol , ySol2]
    label = ['1','2']
    title = 'Many'
    multiPlot('Earth',ys,label,title)

    


if __name__ == '__main__':
    main()