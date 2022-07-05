from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt
import os 
import csv


def multiPlot(NameP,ySol,Labels,Title):
#Ex:ySol = [ySol1 , ySol2]
#Labels and Title must be strings
    N = len(ySol)
    fileName = os.path.join("Data","SolarSystem.csv")
    file = open(fileName,newline='')
    reader = csv.reader(file)
    header = next(reader)
    for row in reader:
        if(row[0] == ''):
            break
        Name = row[0]
        if Name == NameP:
            radius = float(row[2])
    file.close()

    plt.style.use('dark_background')
    fig = plt.figure(figsize=(18,6))
    ax = fig.add_subplot(111,projection='3d')

    for i in range(N):
        ax.plot(ySol[i][:,0],ySol[i][:,1],ySol[i][:,2],label=Labels[i])
        ax.plot([ySol[i][0,0]],[ySol[i][0,1]],[ySol[i][0,2]],'wo')

    _u,_v = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
    _x = radius * np.cos(_u)*np.sin(_v)
    _y = radius * np.sin(_u)*np.sin(_v)
    _z = radius * np.cos(_v)
    ax.plot_surface(_x,_y,_z,cmap='Blues')

    l = radius*2
    x,y,z = [[0,0,0],[0,0,0],[0,0,0]]
    u,v,w = [[l,0,0],[0,l,0],[0,0,l]]
    ax.quiver(x,y,z,u,v,w,color='k')

    maxVal =np.max(np.abs(ySol))

    ax.set_xlim([-maxVal,maxVal])
    ax.set_ylim([-maxVal,maxVal])
    ax.set_zlim([-maxVal,maxVal])

    ax.set_xlabel(['X (km)'])
    ax.set_ylabel(['Y (km)'])
    ax.set_zlabel(['X (km)'])

    ax.set_title(Title)
    plt.legend()
    plt.show()