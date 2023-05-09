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
        line1 = Line([0,0,0], [1,0,0], color=BLUE)
        line2 = Line([0,0,0], [0,-0.5,-0.866], color=RED)
        line3 = Line([0,0,0], [0, 1, 0], color=GREEN)

        self.counter = 0
        self.dt = 0
        self.phi = 0
        self.psi = 0
        self.theta = 0
        self.baseVec1 = line1.get_unit_vector()
        self.baseVec2 = line2.get_unit_vector()
        self.baseVec3 = line3.get_unit_vector()
        self.phiTarget = 0
        self.psiTarget = 118.31816295961119
        self.thetaTarget = 26.931041831092045
        def updater1(mob,dt):
            deg2rad = np.pi/180
            rad2deg = 180/np.pi
            phiGoal = abs(self.phiTarget)
            thetaGoal = abs(self.thetaTarget)
            psiGoal = abs(self.psiTarget)
            self.dt += dt
            if self.dt >= 0.1:
                self.counter += 1
                if self.counter >= 360:
                    self.counter = 0
                    self.phi = 0
                    self.psi = 0
                    self.theta = 0
                self.dt = 0
                if (self.counter*deg2rad) <= (psiGoal*deg2rad):
                    if self.psiTarget >= 0:
                        self.psi += (1*deg2rad)
                    else:
                        self.psi -= (1*deg2rad)
                if ((self.counter*deg2rad) >= (psiGoal*deg2rad)) and((self.counter*deg2rad) <= (psiGoal*deg2rad + thetaGoal*deg2rad)):
                    if self.thetaTarget >= 0:
                        self.theta += (1*deg2rad)
                    else:
                        self.theta -= (1*deg2rad)
                if ((self.counter*deg2rad) >= (phiGoal*deg2rad)) and ((self.counter*deg2rad) > (psiGoal*deg2rad + thetaGoal*deg2rad)) and((self.counter*deg2rad) <= (psiGoal*deg2rad + phiGoal*deg2rad + thetaGoal*deg2rad)):
                    if self.phiTarget >= 0:
                        self.phi += (1*deg2rad)
                    else:
                        self.phi -= (1*deg2rad)
                
                Rx = [1, 0, 0,
                    0, m.cos(self.phi), -m.sin(self.phi),
                    0, m.sin(self.phi), m.cos(self.phi)]
                Ry = [m.cos(self.theta), 0, m.sin(self.theta),
                    0, 1, 0,
                    -m.sin(self.theta), 0, m.cos(self.theta)]
                Rz = [m.cos(self.psi), -m.sin(self.psi), 0,
                    m.sin(self.psi), m.cos(self.psi), 0,
                    0, 0, 1]
                Rx = np.array(Rx)
                Rx = Rx.reshape((3,3))
                Ry = np.array(Ry)
                Ry = Ry.reshape((3,3))
                Rz = np.array(Rz)
                Rz = Rz.reshape((3,3))
                DCM = np.matmul(Ry,Rx)
                DCM = np.matmul(Rz,DCM)
                r_prime = np.matmul(DCM,self.baseVec1)
                line = Line([0,0,0], r_prime, color=BLUE)
                print(r_prime)
                mob.become(line)

        line1.add_updater(updater1)
        self.add(ax, line1, line2,line3)
        self.wait(60)
