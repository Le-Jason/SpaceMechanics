import sys
sys.path.append('C:\\Users\\spong\\Documents\\Code\\projects\\manimrepo\\manim')
from manimlib import *
import numpy as np
sys.path.append('C:/Users/spong/Documents/Code/projects/spacemechanics/tools/data')
sys.path.append('C:/Users/spong/Documents/Code/projects/spacemechanics/tools')
from planetaryData import earth
from spaceCraftState import spaceCraftState
from spaceCraftState import perturbations
import math as m
from numpy.linalg import inv

class Vectors(Scene):
    def construct(self):
        frame=self.camera.frame
        frame.set_euler_angles(
            theta=50 * DEGREES,
            phi=60 * DEGREES,
        )
        frame.scale(0.4)
        ax = ThreeDAxes()
        scaleFactor = 1
        vector1 = [898.357/scaleFactor,-330.16/scaleFactor,190.49/scaleFactor]
        vector2 = [902.3/scaleFactor,-75.827/scaleFactor,43.425/scaleFactor]
        # vector3a = [vector1[0]-vector2[0],0,0]
        # vector3b = [0,vector1[1]-vector2[1],0]
        # vector3c = [0,0,vector1[2]-vector2[2]]
        vector3 = [vector1[0]-vector2[0],vector1[1]-vector2[1],vector1[2]-vector2[2]]
        print(vector3)
        deg2rad = np.pi/180
        phi = -30*deg2rad
        theta = 0
        psi = 0
        Rx = [1, 0, 0,
            0, m.cos(phi), -m.sin(phi),
            0, m.sin(phi), m.cos(phi)]
        Ry = [m.cos(theta), 0, m.sin(theta),
            0, 1, 0,
            -m.sin(theta), 0, m.cos(theta)]
        Rz = [m.cos(psi), -m.sin(psi), 0,
            m.sin(psi), m.cos(psi), 0,
            0, 0, 1]
        Rx = np.array(Rx)
        Rx = Rx.reshape((3,3))
        Ry = np.array(Ry)
        Ry = Ry.reshape((3,3))
        Rz = np.array(Rz)
        Rz = Rz.reshape((3,3))
        DCM = np.matmul(Ry,Rz)
        DCM = np.matmul(Rx,DCM)
        r_prime = np.matmul(Rx,vector3)
        print(vector3)
        print(r_prime)
        vector3b = [r_prime[0],0,0]
        vector3a = [0,r_prime[1],0]
        vector3c = [0,0,r_prime[2]]




        line1 = Line([0,0,0], vector1, color=BLUE)
        line2 = Line([0,0,0], vector2, color=RED)
        line3 = Line([0,0,0], r_prime, color=GREEN)
        line3a = Line([0,0,0], vector3a, color=YELLOW)
        line3b = Line([0,0,0], vector3b, color=YELLOW)
        line3c = Line([0,0,0], vector3c, color=YELLOW)
        
        self.counter1 = 0
        self.dt1 = 0
        self.phi1 = 0
        self.psi1 = 0
        self.theta1 = 0
        self.counter2 = 0
        self.dt2 = 0
        self.phi2 = 0
        self.psi2 = 0
        self.theta2 = 0
        self.counter3 = 0
        self.dt3 = 0
        self.phi3 = 0
        self.psi3 = 0
        self.theta3 = 0
        self.baseVec1 = vector1
        self.baseVec2 = vector2
        self.baseVec3a = vector3a
        self.baseVec3b = vector3b
        self.baseVec3c = vector3c
        self.phiTarget = 180
        self.psiTarget = 30
        self.thetaTarget = -10
        def updater1(mob,dt):
            deg2rad = np.pi/180
            rad2deg = 180/np.pi
            phiGoal = abs(self.phiTarget)
            thetaGoal = abs(self.thetaTarget)
            psiGoal = abs(self.psiTarget)
            self.dt1 += dt
            if self.dt1 >= 0.1:
                self.counter1 += 1
                if self.counter1 >= 360:
                    self.counter1 = 0
                    self.phi1 = 0
                    self.psi1 = 0
                    self.theta1 = 0
                self.dt1 = 0
                if (self.counter1*deg2rad) <= (psiGoal*deg2rad):
                    if self.psiTarget >= 0:
                        self.psi1 += (1*deg2rad)
                    else:
                        self.psi1 -= (1*deg2rad)
                if ((self.counter1*deg2rad) >= (psiGoal*deg2rad)) and((self.counter1*deg2rad) <= (psiGoal*deg2rad + thetaGoal*deg2rad)):
                    if self.thetaTarget >= 0:
                        self.theta1 += (1*deg2rad)
                    else:
                        self.theta1 -= (1*deg2rad)
                if ((self.counter1*deg2rad) > (psiGoal*deg2rad + thetaGoal*deg2rad)) and((self.counter1*deg2rad) <= (psiGoal*deg2rad + phiGoal*deg2rad + thetaGoal*deg2rad)):
                    if self.phiTarget >= 0:
                        self.phi1 += (1*deg2rad)
                    else:
                        self.phi1 -= (1*deg2rad)
                
                Rx = [1, 0, 0,
                    0, m.cos(self.psi1), -m.sin(self.psi1),
                    0, m.sin(self.psi1), m.cos(self.psi1)]
                Ry = [m.cos(self.phi1), 0, m.sin(self.phi1),
                    0, 1, 0,
                    -m.sin(self.phi1), 0, m.cos(self.phi1)]
                Rz = [m.cos(self.theta1), -m.sin(self.theta1), 0,
                    m.sin(self.theta1), m.cos(self.theta1), 0,
                    0, 0, 1]
                Rx = np.array(Rx)
                Rx = Rx.reshape((3,3))
                Ry = np.array(Ry)
                Ry = Ry.reshape((3,3))
                Rz = np.array(Rz)
                Rz = Rz.reshape((3,3))
                DCM = np.matmul(Rz,Ry)
                DCM = np.matmul(Rx,DCM)
                self.r_prime = np.matmul(DCM,self.baseVec3a)
                line = Line([0,0,0], self.r_prime, color=YELLOW)
                mob.become(line)
        def updater2(mob,dt):
            deg2rad = np.pi/180
            rad2deg = 180/np.pi
            phiGoal = abs(self.phiTarget)
            thetaGoal = abs(self.thetaTarget)
            psiGoal = abs(self.psiTarget)
            self.dt2 += dt
            if self.dt2 >= 0.1:
                self.counter2 += 1
                if self.counter2 >= 360:
                    self.counter2 = 0
                    self.phi2 = 0
                    self.psi2 = 0
                    self.theta2 = 0
                self.dt2 = 0
                if (self.counter2*deg2rad) <= (psiGoal*deg2rad):
                    if self.psiTarget >= 0:
                        self.psi2 += (1*deg2rad)
                    else:
                        self.psi2 -= (1*deg2rad)
                if ((self.counter2*deg2rad) >= (psiGoal*deg2rad)) and((self.counter2*deg2rad) <= (psiGoal*deg2rad + thetaGoal*deg2rad)):
                    if self.thetaTarget >= 0:
                        self.theta2 += (1*deg2rad)
                    else:
                        self.theta2 -= (1*deg2rad)
                if ((self.counter2*deg2rad) > (psiGoal*deg2rad + thetaGoal*deg2rad)) and((self.counter2*deg2rad) <= (psiGoal*deg2rad + phiGoal*deg2rad + thetaGoal*deg2rad)):
                    if self.phiTarget >= 0:
                        self.phi2 += (1*deg2rad)
                    else:
                        self.phi2 -= (1*deg2rad)
                
                Rx = [1, 0, 0,
                    0, m.cos(self.psi2), -m.sin(self.psi2),
                    0, m.sin(self.psi2), m.cos(self.psi2)]
                Ry = [m.cos(self.phi2), 0, m.sin(self.phi2),
                    0, 1, 0,
                    -m.sin(self.phi2), 0, m.cos(self.phi2)]
                Rz = [m.cos(self.theta2), -m.sin(self.theta2), 0,
                    m.sin(self.theta2), m.cos(self.theta2), 0,
                    0, 0, 1]
                Rx = np.array(Rx)
                Rx = Rx.reshape((3,3))
                Ry = np.array(Ry)
                Ry = Ry.reshape((3,3))
                Rz = np.array(Rz)
                Rz = Rz.reshape((3,3))
                DCM = np.matmul(Rz,Ry)
                DCM = np.matmul(Rx,DCM)
                r_prime = np.matmul(DCM,self.baseVec3a)
                r_prime2 = np.matmul(DCM,self.baseVec3b)
                adder = [r_prime[0]+r_prime2[0],r_prime[1]+r_prime2[1],r_prime[2]+r_prime2[2]]
                line = Line(r_prime, adder, color=YELLOW)
                mob.become(line)
        def updater3(mob,dt):
            deg2rad = np.pi/180
            rad2deg = 180/np.pi
            phiGoal = abs(self.phiTarget)
            thetaGoal = abs(self.thetaTarget)
            psiGoal = abs(self.psiTarget)
            self.dt3 += dt
            if self.dt3 >= 0.1:
                self.counter3 += 1
                if self.counter3 >= 360:
                    self.counter3 = 0
                    self.phi3 = 0
                    self.psi3 = 0
                    self.theta3 = 0
                self.dt3 = 0
                if (self.counter3*deg2rad) <= (psiGoal*deg2rad):
                    if self.psiTarget >= 0:
                        self.psi3 += (1*deg2rad)
                    else:
                        self.psi3 -= (1*deg2rad)
                if ((self.counter3*deg2rad) >= (psiGoal*deg2rad)) and((self.counter3*deg2rad) <= (psiGoal*deg2rad + thetaGoal*deg2rad)):
                    if self.thetaTarget >= 0:
                        self.theta3 += (1*deg2rad)
                    else:
                        self.theta3 -= (1*deg2rad)
                if ((self.counter3*deg2rad) > (psiGoal*deg2rad + thetaGoal*deg2rad)) and((self.counter3*deg2rad) <= (psiGoal*deg2rad + phiGoal*deg2rad + thetaGoal*deg2rad)):
                    if self.phiTarget >= 0:
                        self.phi3 += (1*deg2rad)
                    else:
                        self.phi3 -= (1*deg2rad)
                
                Rx = [1, 0, 0,
                    0, m.cos(self.psi3), -m.sin(self.psi3),
                    0, m.sin(self.psi3), m.cos(self.psi3)]
                Ry = [m.cos(self.phi3), 0, m.sin(self.phi3),
                    0, 1, 0,
                    -m.sin(self.phi3), 0, m.cos(self.phi3)]
                Rz = [m.cos(self.theta3), -m.sin(self.theta3), 0,
                    m.sin(self.theta3), m.cos(self.theta3), 0,
                    0, 0, 1]
                Rx = np.array(Rx)
                Rx = Rx.reshape((3,3))
                Ry = np.array(Ry)
                Ry = Ry.reshape((3,3))
                Rz = np.array(Rz)
                Rz = Rz.reshape((3,3))
                DCM = np.matmul(Rz,Ry)
                DCM = np.matmul(Rx,DCM)
                r_prime = np.matmul(DCM,self.baseVec3a)
                r_prime2 = np.matmul(DCM,self.baseVec3b)
                r_prime3 = np.matmul(DCM,self.baseVec3c)
                adder = [r_prime[0]+r_prime2[0],r_prime[1]+r_prime2[1],r_prime[2]+r_prime2[2]]
                adder2 = [r_prime[0]+r_prime2[0]+r_prime3[0],r_prime[1]+r_prime2[1]+r_prime3[1],r_prime[2]+r_prime2[2]+r_prime3[2]]
                print(self.counter3)
                print(self.phi3*180/np.pi)
                line = Line(adder, adder2, color=YELLOW)
                mob.become(line)


        # line3a.add_updater(updater1)
        # line3b.add_updater(updater2)
        # line3c.add_updater(updater3)
        self.add(ax, line1, line2,line3a,line3b,line3c)
        self.wait(60)
