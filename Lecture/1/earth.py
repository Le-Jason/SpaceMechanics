import PIL
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import sys
plt.rcParams['animation.ffmpeg_path'] = 'C:\\Users\\spong\\Downloads\\ffmpeg-6.0-full_build\\ffmpeg-6.0-full_build\\bin\\ffmpeg.exe'


sys.path.append('C:/Users/spong/Documents/Code/projects/spacemechanics/tools/data')
sys.path.append('C:/Users/spong/Documents/Code/projects/spacemechanics/tools')
from planetaryData import earth
from spaceCraftState import spaceCraftState
from spaceCraftState import perturbations

# load bluemarble with PIL
bm = PIL.Image.open('earthmap.jpg')
# it's big, so I'll rescale it, convert to array, and divide by 256 to get RGB values that matplotlib accept 
bm = np.array(bm.resize([int(d/5) for d in bm.size]))/256.

# coordinates of the image - don't know if this is entirely accurate, but probably close
lons = np.linspace(-180, 180, bm.shape[1]) * np.pi/180 
lats = np.linspace(-90, 90, bm.shape[0])[::-1] * np.pi/180 

# repeat code from one of the examples linked to in the question, except for specifying facecolors:
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d',facecolor='black')
plt.style.use('dark_background')

state = [414+6371+1000+10000,0.5006189,0,105.6372,234.1955,0]
perts=perturbations()
perts['J2'] + True
Orbit1 = spaceCraftState('Earth',kepler=True,stateVec=state,perturbation=perts)
t0 = 0
y0 = [0,0,0,0,0,0]
dt = 100
tf = 100*5000
Labels = ['GPS BIIR-2']
Title = 'GPS'
tSol,ySol = Orbit1.propagateOrbit(t0,y0,dt,tf)
N = 1
# for i in range(N):
#     ax.plot(ySol[:,0],ySol[:,1],ySol[:,2],label=Labels[i],linewidth=1,color='purple')
#     ax.plot([ySol[0,0]],[ySol[0,1]],[ySol[0,2]],'wo')

maxVal =np.max(np.abs(ySol))

print(ySol)


x = np.outer(np.cos(lons), np.cos(lats)).T
y = np.outer(np.sin(lons), np.cos(lats)).T
z = np.outer(np.ones(np.size(lons)), np.sin(lats)).T
ax.plot_surface(x*earth['radius'], y*earth['radius'], z*earth['radius'], rstride=4, cstride=4, facecolors = bm)

xlist = []
y1list = []
y2list = []
y3list = []

metaData = dict(title='Movie',artist='FUCKYOU')
writer = FFMpegWriter(fps=30,metadata=metaData)

with writer.saving(fig,"fig.mp4",100):
    # for i in range(N):
    #     ax.plot(ySol[:,0],ySol[:,1],ySol[:,2],label=Labels[i],linewidth=1,color='purple')
    #     ax.plot([ySol[0,0]],[ySol[0,1]],[ySol[0,2]],'wo')
    # x = np.outer(np.cos(lons), np.cos(lats)).T
    # y = np.outer(np.sin(lons), np.cos(lats)).T
    # z = np.outer(np.ones(np.size(lons)), np.sin(lats)).T
    # ax.plot_surface(x*earth['radius'], y*earth['radius'], z*earth['radius'], rstride=4, cstride=4, facecolors = bm)
    # for i in range(len(tSol)):
    for i in range(len(tSol[0:300])):
        ax.plot_surface(x*earth['radius'], y*earth['radius'], z*earth['radius'], rstride=4, cstride=4, facecolors = bm)
        xlist.append(tSol[i])
        y1list.append(ySol[i,0])
        y2list.append(ySol[i,1])
        y3list.append(ySol[i,2])
        ax.set_xlim([-maxVal,maxVal])
        ax.set_ylim([-maxVal,maxVal])
        ax.set_zlim([-maxVal,maxVal])
        ax.plot(y1list,y2list,y3list,color='black')
        ax.scatter(ySol[i,0],ySol[i,1],ySol[i,2],color='red', linewidth=1)
        writer.grab_frame()
        plt.cla()

# plt.show()
