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

class Vectors(Scene):
    def construct(self):
        frame=self.camera.frame
        frame.set_euler_angles(
            theta=50 * DEGREES,
            phi=60 * DEGREES,
        )
        frame.scale(0.4)
        ax = ThreeDAxes()
        line1 = Line([0,0,0], [0,1,0], color=RED)
        line2 = Line([0,0,0], [1,0,0], color=RED)
        line3 = Line([0,0,0], [0,0,1], color=RED)
        self.counter = 0
        self.dt = 0
        self.phi = 0
        self.psi = 0
        self.theta = 0
        self.baseVec1 = line1.get_unit_vector()
        self.baseVec2 = line2.get_unit_vector()
        self.baseVec3 = line3.get_unit_vector()
        def updater1(mob,dt):
            deg2rad = np.pi/180
            rad2deg = 180/np.pi
            phiGoal = 30
            thetaGoal = 180
            psiGoal = 30
            self.dt += dt
            if self.dt >= 0.1:
                self.counter += 1
                if self.counter >= 360:
                    self.counter = 0
                    self.phi = 0
                    self.psi = 0
                    self.theta = 0
                self.dt = 0
                if (self.counter*deg2rad) <= (phiGoal*deg2rad):
                    self.phi += (1*deg2rad)
                if ((self.counter*deg2rad) >= (phiGoal*deg2rad)) and((self.counter*deg2rad) <= (psiGoal*deg2rad + phiGoal*deg2rad)):
                    self.psi += (1*deg2rad)
                if ((self.counter*deg2rad) >= (thetaGoal*deg2rad)) and ((self.counter*deg2rad) > (psiGoal*deg2rad + phiGoal*deg2rad)) and((self.counter*deg2rad) <= (psiGoal*deg2rad + phiGoal*deg2rad + thetaGoal*deg2rad)):
                    self.theta += (1*deg2rad)
                # if (self.counter*deg2rad) <= (psiGoal*deg2rad):
                #     self.psi += (1*deg2rad)
                # if ((self.counter*deg2rad) >= (psiGoal*deg2rad)) and((self.counter*deg2rad) <= (phiGoal*deg2rad + psiGoal*deg2rad)):
                #     self.phi += (1*deg2rad)
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
                DCM = np.matmul(Rz,Ry)
                DCM = np.matmul(Rx,DCM)
                r_prime = np.matmul(DCM,self.baseVec1)
                line = Line([0,0,0], r_prime, color=RED)
                mob.become(line)
            
        def updater2(mob,dt):
            deg2rad = np.pi/180
            rad2deg = 180/np.pi
            phiGoal = 30
            thetaGoal = 180
            psiGoal = 30
            self.dt += dt
            if self.dt >= 0.1:
                self.counter += 1
                if self.counter >= 360:
                    self.counter = 0
                    self.phi = 0
                    self.psi = 0
                    self.theta = 0
                self.dt = 0
                if (self.counter*deg2rad) <= (phiGoal*deg2rad):
                    self.phi += (1*deg2rad)
                if ((self.counter*deg2rad) >= (phiGoal*deg2rad)) and((self.counter*deg2rad) <= (psiGoal*deg2rad + phiGoal*deg2rad)):
                    self.psi += (1*deg2rad)
                if ((self.counter*deg2rad) >= (thetaGoal*deg2rad)) and ((self.counter*deg2rad) > (psiGoal*deg2rad + phiGoal*deg2rad)) and((self.counter*deg2rad) <= (psiGoal*deg2rad + phiGoal*deg2rad + thetaGoal*deg2rad)):
                    self.theta += (1*deg2rad)
                # if (self.counter*deg2rad) <= (psiGoal*deg2rad):
                #     self.psi += (1*deg2rad)
                # if ((self.counter*deg2rad) >= (psiGoal*deg2rad)) and((self.counter*deg2rad) <= (phiGoal*deg2rad + psiGoal*deg2rad)):
                #     self.phi += (1*deg2rad)
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
                DCM = np.matmul(Rz,Ry)
                DCM = np.matmul(Rx,DCM)
                r_prime = np.matmul(DCM,self.baseVec2)
                line = Line([0,0,0], r_prime, color=RED)
                mob.become(line)
        def updater3(mob,dt):
            deg2rad = np.pi/180
            rad2deg = 180/np.pi
            phiGoal = 30
            thetaGoal = 180
            psiGoal = 30
            self.dt += dt
            if self.dt >= 0.1:
                self.counter += 1
                if self.counter >= 360:
                    self.counter = 0
                    self.phi = 0
                    self.psi = 0
                    self.theta = 0
                self.dt = 0
                if (self.counter*deg2rad) <= (phiGoal*deg2rad):
                    self.phi += (1*deg2rad)
                if ((self.counter*deg2rad) >= (phiGoal*deg2rad)) and((self.counter*deg2rad) <= (psiGoal*deg2rad + phiGoal*deg2rad)):
                    self.psi += (1*deg2rad)
                if ((self.counter*deg2rad) >= (thetaGoal*deg2rad)) and ((self.counter*deg2rad) > (psiGoal*deg2rad + phiGoal*deg2rad)) and((self.counter*deg2rad) <= (psiGoal*deg2rad + phiGoal*deg2rad + thetaGoal*deg2rad)):
                    self.theta += (1*deg2rad)
                # if (self.counter*deg2rad) <= (psiGoal*deg2rad):
                #     self.psi += (1*deg2rad)
                # if ((self.counter*deg2rad) >= (psiGoal*deg2rad)) and((self.counter*deg2rad) <= (phiGoal*deg2rad + psiGoal*deg2rad)):
                #     self.phi += (1*deg2rad)
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
                DCM = np.matmul(Rz,Ry)
                DCM = np.matmul(Rx,DCM)
                r_prime = np.matmul(DCM,self.baseVec3)
                line = Line([0,0,0], r_prime, color=RED)
                mob.become(line)

        line1.add_updater(updater1)
        line2.add_updater(updater2)
        line3.add_updater(updater3)
        self.add(ax, line1, line2, line3)

        
        self.wait(60)
