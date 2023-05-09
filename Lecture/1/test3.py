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
        line1 = Line([0,0,0], [-0.423,0.785,-0.453], color=RED)
        line2 = Line([0,0,0], [0,-0.5,-0.866], color=RED)
        line3 = Line([0,0,0], np.cross([-0.423,0.785,-0.453],[0,-0.5,-0.866]), color=RED)

        firstAngle = 90-(m.acos(abs(np.dot([-0.423,0.785,-0.453],[0,0,1])/m.sqrt(0.423**2 + 0.785**2 + 0.453**2)))*180/np.pi)
        print(firstAngle)

        # Define the vector
        v = np.array([-0.423, 0.785, -0.453])

        # Define the xy plane normal vector
        n = np.array([0, 0, 1])

        # Find the angle between the vector and the xy plane
        theta = np.arccos(np.dot(v, n)/(np.linalg.norm(v)*np.linalg.norm(n)))
        

        # Define the rotation matrix
        R = np.array([[np.cos(theta), np.sin(theta), 0],
                    [-np.sin(theta), np.cos(theta), 0],
                    [0, 0, 1]])

        # Rotate the vector
        v_rotated = np.dot(R, v)

        # Project the rotated vector onto the xy plane
        v_projected = np.array([v_rotated[0], v_rotated[1], 0])

        # Normalize the projected vector
        v_normalized = v_projected/np.linalg.norm(v_projected)

        # Check that the projected vector is in the xy plane
        assert np.isclose(np.dot(v_normalized, n), 0)

        # Print the resulting vector in the xy plane
        print(v_normalized* np.linalg.norm(v_projected))

        self.counter = 0
        self.dt = 0
        self.phi = 0
        self.psi = 0
        self.theta = 0
        self.baseVec1 = line1.get_unit_vector()
        self.baseVec2 = line2.get_unit_vector()
        self.baseVec3 = line3.get_unit_vector()
        self.phiTarget = 0
        self.psiTarget = theta*180/np.pi
        self.thetaTarget = 0
        def updater1(mob,dt):
            deg2rad = np.pi/180
            rad2deg = 180/np.pi
            phiGoal = self.phiTarget
            thetaGoal = self.thetaTarget
            psiGoal = self.psiTarget
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
                    self.psi += (1*deg2rad)

                # Rx = [1, 0, 0,
                #     0, m.cos(self.phi), -m.sin(self.phi),
                #     0, m.sin(self.phi), m.cos(self.phi)]
                # Ry = [m.cos(self.theta), 0, m.sin(self.theta),
                #     0, 1, 0,
                #     -m.sin(self.theta), 0, m.cos(self.theta)]
                Rz = [m.cos(self.psi), -m.sin(self.psi), 0,
                    m.sin(self.psi), m.cos(self.psi), 0,
                    0, 0, 1]
                # Rx = np.array(Rx)
                # Rx = Rx.reshape((3,3))
                # Ry = np.array(Ry)
                # Ry = Ry.reshape((3,3))
                Rz = np.array(Rz)
                Rz = Rz.reshape((3,3))
                # DCM = np.matmul(Ry,Rx)
                DCM = Rz
                print(DCM)
                r_prime = np.matmul(DCM,self.baseVec1)
                line = Line([0,0,0], r_prime, color=RED)
                mob.become(line)

        line1.add_updater(updater1)
        self.add(ax, line1, line2, line3)

        
        self.wait(60)
