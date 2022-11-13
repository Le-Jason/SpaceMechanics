from vpython import *
import numpy as np
import time
import math
import matplotlib.pyplot as plt
import matplotlib.animation

def eulerAngleRotation(phi,theta,psi,vector):
    R_x = np.array([[1, 0, 0],
                    [0, math.cos(phi), -math.sin(phi) ],
                    [0, math.sin(phi), math.cos(phi)  ]
                    ])
 
    R_y = np.array([[math.cos(theta), 0, math.sin(theta)  ],
                    [0, 1, 0],
                    [-math.sin(theta), 0, math.cos(theta)  ]
                    ])
 
    R_z = np.array([[math.cos(psi), -math.sin(psi), 0],
                    [math.sin(psi), math.cos(psi), 0],
                    [0, 0, 1]
                    ])
 
    R = np.dot(R_z, np.dot( R_y, R_x ))
    return np.matmul(vector,R)

def Ysolution(p,n,coeff,t):
    x0 = p[0]
    z0 = p[2]
    a = n[0]
    c = n[2]
    d = coeff
 
    x = np.zeros((len(t),1))
    y = np.zeros((len(t),1))
    z = np.zeros((len(t),1))

    x = []
    y = []
    z = []

    CNT = 0
    for i in t:
        try:
            y.append(math.sqrt((((((-a/c)*(i-x0)) + z0)**2)/d)-(i**2)))
            x.append(i)
            z.append(((-a/c)*(i-x0)) + z0)
        except:
            pass
        CNT += 1
    return x,y,z

def MainPolySol(a,c,d,x0,z0):
    Part1 = (a**2)/((c**2)*d)
    Part2 = -(2*(a**2)*x0)/((c**2)*d)
    Part3 = ((a**2) * (x0**2))/((c**2)*d)
    Part4 = -(2*a*z0)/(c*d)
    Part5 = (2*a*x0*z0)/(c*d)
    Part6 = (z0**2)/d
    Part7 = -1
    g = Part1 + Part7
    h = Part2 + Part4
    f = Part3 + Part5 + Part6
    tOne = (-h + math.sqrt((h**2) -(4*g*f)))/(2*g)
    tTwo = (-h - math.sqrt((h**2) -(4*g*f)))/(2*g)

    Part1 = (2*(a**2))/((c**2)*d)
    Part2 = (-2*(a**2)*x0)/((c**2)*d)
    Part3 = (-2*a*z0)/(c*d)
    Part4 = -2

    Part1 = a*((a*x0)+(c*z0))/((a**2)-((c**2)*d))
    point = Part1

    print("This is plus:" + str(tOne))
    print("This is neg:" + str(tTwo))


    Part1 = (2*(a**2))/((c**2)*d) -2
    val = Part1
    N = 10000
    if val < 0:
        t = np.linspace(tOne,tTwo,num=N)
    else:
        arr = [tOne,tTwo]
        minT = min(arr)
        maxT = max(arr)
   
        tG = np.linspace(-minT*4,-minT,num=5000)
        bG = np.linspace(maxT,maxT*4,num=5000)
        t = np.concatenate((tG,bG),axis=None)
    
    return t

plt.style.use('dark_background')
fig = plt.figure(figsize=(18,6))
ax = fig.add_subplot(111,projection='3d')
ax.set_xlabel(["X (km)"])
ax.set_ylabel(["Y (km)"])
ax.set_zlabel(["Z (km)"])
ax.set_title('Title')

deg2rad = np.pi/180
n = [0,0,1]

p = [-1,0,-1]
coeff = 4
t = np.linspace(-1.5,1.5,num=500)
x = 0
y = 0
z = 0
line1, = ax.plot(x,y)

scene.range=5
scene.forward=vector(-3,-3,-3)
 
scene.width=600
scene.height=600
 
bottomCone = cone(pos=vector(0,-5,0),axis=vector(0,1,0),size=vector(5,5,5),color=color.white,opacity=0.4)
topCone = cone(pos=vector(0,5,0),axis=vector(0,-1,0),size=vector(5,5,5),color=color.white,opacity=0.4)
plane = box(pos=vector(-1,-1,0),axis=vector(1,0,0),size=vector(20,0.01,20),color=color.blue,opacity=0.5)

deg2rad = np.pi/180

prime = [1,0,0]


def update(i):

    print(i)

    [a,b,c] = eulerAngleRotation(0,0,i*deg2rad,prime)

    plane.axis = vector(a,b,c)
    plane.size=vector(20,0.1,20)

    x0=-1
    z0=-1
    d=4
    if i < 45:
        time.sleep(0.2)
    else:
        time.sleep(0.2)
    nPrime = eulerAngleRotation(0,-i*deg2rad,0,n)
    t = MainPolySol(nPrime[0],nPrime[2],d,x0,z0)
    [x,y,z] = Ysolution(p,nPrime,coeff,t)
    plt.cla()
    yNeg = []

    for i in y:
        yNeg.append(-1*i)

    xSum = 0
    ySum = 0
    zSum = 0
    # for i in range(len(x)):
    #     xSum += x[i]
    #     ySum += y[i]
    #     zSum += z[i]

    plt.plot(x,y,z,'w',label='Trajectory')
    plt.plot(x,yNeg,z,'w',label='Trajectory')
    # print(xSum)
    # print(ySum)
    # print(zSum)
    # plt.plot(xSum/len(x),ySum/len(x),zSum/len(x),'r',label='Trajectory')
    # plt.plot(0,0,0,'r',label='Trajectory')


    plt.tight_layout()
    plt.title('Orbit')

    
    ax.set_xlim([-20,20])
    ax.set_ylim([-20,20])
    ax.set_zlim([-20,20])

    


ani = matplotlib.animation.FuncAnimation(fig, update, frames=90, repeat=True)    
plt.show()
    
