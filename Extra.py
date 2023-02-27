from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt
import math as m
import os 
import csv
from mpl_toolkits.mplot3d import Axes3D

def multiPlot(NameP,ySol,Labels,Title):
#Ex:ySol = [ySol1 , ySol2]
#Labels and Title must be strings
    N = len(ySol)
    fileName = os.path.join("C:/Users/spong/Documents/Code/projects/spacemechanics/Data","SolarSystem.csv")
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

    img = plt.imread('C:/Users/spong/Documents/Code/projects/spacemechanics/Data/bluemarble.jpg')
    theta = np.linspace(0, np.pi, img.shape[0])
    phi = np.linspace(0, 2*np.pi, img.shape[1])
    count = 180
    theta_inds = np.linspace(0, img.shape[0] - 1, count).round().astype(int)
    phi_inds = np.linspace(0, img.shape[1] - 1, count).round().astype(int)
    theta = theta[theta_inds]
    phi = phi[phi_inds]
    img = img[np.ix_(theta_inds, phi_inds)]
    theta,phi = np.meshgrid(theta, phi)

    phi,theta = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
    _x = radius * np.cos(phi)*np.sin(theta)
    _y = radius * np.sin(phi)*np.sin(theta)
    _z = radius * np.cos(theta)
    # ax.plot_surface(_x.T,_y.T,_z.T, facecolors=img/255, cstride=1, rstride=1)
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


def readData(name):
    fileName = os.path.join("C:/Users/spong/Documents/Code/projects/spacemechanics/Data","SolarSystem.csv")
    file = open(fileName,newline='')
    reader = csv.reader(file)
    header = next(reader)
    for row in reader:
        if(row[0] == ''):
            break
        Name = row[0]
        if name == Name:
            mass = float(row[1])
            radius = float(row[2])
            rotation = float(row[3])
            distance = float(row[4])
            eccentricity  = float(row[5])
            J2 = float(row[6])
            file.close()
            return [Name,mass,radius,rotation,distance,eccentricity,J2]
    file.close()
    return 0


def giveOrbitalElements(rSol,vSol,mu):
    n = len(rSol[:,0])
    aElements = []
    eElements = []
    iElements = []
    RAANElements = []
    wElements = []
    vuElements = []

    for i in range(n):
        r = np.array([rSol[i,0],rSol[i,1],rSol[i,2]])
        v = np.array([vSol[i,0],vSol[i,1],vSol[i,2]])
        K = np.array([0,0,1])
        rNorm = m.sqrt(np.dot(r,r))
        vNorm = m.sqrt(np.dot(v,v))
        vRad = np.dot(r,v)/rNorm
        h = np.cross(r,v)
        hNorm = m.sqrt(np.dot(h,h))
        i = m.acos(h[2],hNorm)
        N = np.cross(K,h)
        NNorm = m.sqrt(np.dot(N,N))
        if N[1] >= 0:
            RAAN = m.acos(N[0]/NNorm)
        elif N[1] < 0:
            RAAN = (2*np.pi) - m.acos(N[0]/NNorm)
        e = (1/mu)*(np.cross(v,h) - (mu*(r/rNorm)))
        eNorm = m.sqrt(np.dot(e,e))
        if e[2] >= 0:
            w = m.acos(np.dot(N,e)/(NNorm*eNorm))
        elif e[2] < 0:
            w = (2*np.pi) - m.acos(np.dot(N,e)/(NNorm*eNorm))
        vuPart = np.dot(e/eNorm,r/rNorm)
        tol = 0.00000000002
        if (vuPart > 1) and (vuPart <= 1+tol):
            vuPart = 1
        if vRad >= 0:
            vu = m.acos(vuPart)
        elif vRad < 0:
            vu = (2*np.pi) - m.acos(vuPart)
        rp = ((hNorm**2)/mu*(1/(1+eNorm*m.cos(0))))
        ra = ((hNorm**2)/mu*(1/(1+eNorm*m.cos(m.pi))))
        a = 0.5*(rp + ra)

        aElements.append(a)
        eElements.append(eNorm)
        iElements.append(i)
        RAANElements.append(RAAN)
        wElements.append(w)
        vuElements.append(vu)
    
        return aElements,eElements,iElements,RAANElements,wElements,vuElements