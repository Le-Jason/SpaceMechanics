from scipy.integrate import ode
import numpy as np
import matplotlib.pyplot as plt
import math as m
import os 
import csv

class orbitInitPropagate():
    def __init__(self,planet,kepler='False',stateVec=[]):
        self.name = planet
        self.readData()
        self.G = 0.0000000000000000000667430 #(N*km^2)/kg^2
        self.mu = self.G * self.mass
        self.deg2rad = np.pi/180
        self.rad2deg = 180/np.pi
        self.kepler = kepler

        if kepler:
            self.stateVec = stateVec
            self.a = stateVec[0]  #semi-major axis [km]
            self.e = stateVec[1]  #eccentricity [--]
            self.i = stateVec[2]*self.deg2rad  #inclination [rad]
            self.RAAN = stateVec[3]*self.deg2rad #Right ascension of the ascending node [rad]
            self.w = stateVec[4]*self.deg2rad  #argument of perigee [rad]
            self.vu = stateVec[5]*self.deg2rad #true anomaly [rad]
            self.COE2RV()

    def readData(self):
        fileName = os.path.join("Data","SolarSystem.csv")
        file = open(fileName,newline='')
        reader = csv.reader(file)
        header = next(reader)
        for row in reader:
            if(row[0] == ''):
                break
            Name = row[0]
            if self.name == Name:
                self.mass = float(row[1])
                self.radius = float(row[2])
                self.rotation = float(row[3])
                self.distance = float(row[4])
                self.eccentricity  = float(row[5])
        file.close()

    def dydt(self,t,y):
        [rx,ry,rz,vx,vy,vz] = y
        rVec = np.array([rx,ry,rz])
        rNorm = np.sqrt(((rx**2)+(ry**2)+(rz**2)))
        [ax,ay,az] = (-rVec*self.mu)/(rNorm**3)
        return [vx,vy,vz,ax,ay,az]

    def propagateOrbit(self,t0,y0,dt,tf):
        N = int(np.ceil(tf/dt))
        if self.kepler:
            y0[0:3] = self.rIJK
            y0[3:7] = self.vIJK
        ySol = np.zeros((N,len(y0)))
        tSol = np.zeros((N,t0))
        for i in range(len(y0)):
            ySol[0,i] = y0[i]
        tSol[0] = t0
        i = 1

        solver = ode(self.dydt)
        solver.set_integrator('lsoda')
        solver.set_initial_value(y0,t0)
        while solver.successful() and i < N:
            solver.integrate(solver.t+dt)
            tSol[i] = solver.t
            ySol[i] = solver.y
            i += 1
        
        self.tSol = tSol
        self.ySol = ySol
        return tSol,ySol
    
    def plot(self):
        plt.style.use('dark_background')
        fig = plt.figure(figsize=(18,6))
        ax = fig.add_subplot(111,projection='3d')

        ax.plot(self.ySol[:,0],self.ySol[:,1],self.ySol[:,2],'w',label='Trajectory')
        
        ax.plot([self.ySol[0,0]],[self.ySol[0,1]],[self.ySol[0,2]],'wo',label='Inital Position')

        _u,_v = np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
        _x = self.radius * np.cos(_u)*np.sin(_v)
        _y = self.radius * np.sin(_u)*np.sin(_v)
        _z = self.radius * np.cos(_v)
        ax.plot_surface(_x,_y,_z,cmap='Blues')

        l = self.radius*2
        x,y,z = [[0,0,0],[0,0,0],[0,0,0]]
        u,v,w = [[l,0,0],[0,l,0],[0,0,l]]
        ax.quiver(x,y,z,u,v,w,color='k')

        maxVal =np.max(np.abs(self.ySol))

        ax.set_xlim([-maxVal,maxVal])
        ax.set_ylim([-maxVal,maxVal])
        ax.set_zlim([-maxVal,maxVal])

        ax.set_xlabel(["X (km)"])
        ax.set_ylabel(["Y (km)"])
        ax.set_zlabel(["Z (km)"])

        ax.set_title('Title')
        plt.legend()
        plt.show()
    
    def COE2RV(self):
        p = self.a*(1 - self.e**2)
        rPQW = np.array([(p*m.cos(self.vu))/(1 + (self.e*m.cos(self.vu))),
                (p*m.sin(self.vu))/(1 + (self.e*m.cos(self.vu))),
                 0])
        rPQW = rPQW.reshape((3,1))

        vPQW = np.array([-m.sqrt(self.mu/p)*m.sin(self.vu),
                m.sqrt(self.mu/p)*(self.e + m.cos(self.vu)),
                0])
        vPQW = vPQW.reshape((3,1))
        
        T_PQW_IJK = [m.cos(self.RAAN)*m.cos(self.w) - m.sin(self.RAAN)*m.sin(self.w)*m.cos(self.i), -m.cos(self.RAAN)*m.sin(self.w) - m.sin(self.RAAN)*m.cos(self.w)*m.cos(self.i), m.sin(self.RAAN)*m.sin(self.i),
                     m.sin(self.RAAN)*m.cos(self.w) + m.cos(self.RAAN)*m.sin(self.w)*m.cos(self.i), -m.sin(self.RAAN)*m.sin(self.w) + m.cos(self.RAAN)*m.cos(self.w)*m.cos(self.i), -m.cos(self.RAAN)*m.sin(self.i),
                     m.sin(self.w)*m.sin(self.i), m.cos(self.w)*m.sin(self.i), m.cos(self.i)]
        T_PQW_IJK = np.array(T_PQW_IJK)
        T_PQW_IJK = T_PQW_IJK.reshape((3,3))


        self.rIJK = np.matmul(T_PQW_IJK,rPQW)
        self.rIJK = self.rIJK.ravel()
        self.vIJK = np.matmul(T_PQW_IJK,vPQW)
        self.vIJK = self.vIJK.ravel()