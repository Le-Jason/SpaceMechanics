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


# Define the vector
v = np.array([-0.423, 0.785, -0.453])

# Define the xy plane normal vector
n = np.array([0, 0, 1])

# Find the angle between thcose vector and the xy plane
theta = abs(np.arccos(np.dot(v, n)/(np.linalg.norm(v)*np.linalg.norm(n))))*180/np.pi
theta = 90-theta
print(theta)
# Define the xz plane normal vector
v = np.array([-0.423, 0.785, 0])
n = np.array([0, 1, 0])

# Find the angle between the vector and the xy plane
psi = abs(np.arccos(np.dot(v, n)/(np.linalg.norm(v)*np.linalg.norm(n))))*180/np.pi
psi = psi + 90
psi = psi*-1
print(psi)
phi = 0

u = np.array([1, 0, 0])
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
DCM = np.matmul(Ry,Rx)
DCM = np.matmul(Rz,DCM)
r_prime = DCM @ u
print(r_prime)
