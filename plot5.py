import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

nParticles = 1000
nSteps = 1000

LoadnParticles = nParticles

coords = []
for particle in range(0,LoadnParticles):
    coords.append(np.zeros([nSteps,3]))
#coord[particleid][step#][xyz=012]


for step in range(0,nSteps):
    datafile = pd.read_csv("data/step%05d.csv"%step,sep=",")

    for particle in range(0,LoadnParticles):
        coords[particle][step][0] = datafile["xPos"][particle]
        coords[particle][step][1] = datafile["yPos"][particle]
        coords[particle][step][2] = datafile["zPos"][particle]



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

steptoplotto = 25
for step in range(0,steptoplotto):
    for particle in range(0,LoadnParticles):
        xPos = coords[particle][step][0]
        yPos = coords[particle][step][1]
        zPos = coords[particle][step][2]
        ax.scatter(xPos, yPos, zPos, zdir='z', s=20, c='k', depthshade=True)
        ax.set_xlim(0, 1024)
        ax.set_ylim(0, 1024)
        ax.set_zlim(0, 1024)

        ax.set_xlabel('X ')
        ax.set_ylabel('Y ')
        ax.set_zlabel('Z ')
    plt.savefig("plots/step%05d.png"%step)
    plt.cla()


'''#Plotting particles in time
for particle in range(0,LoadnParticles):
    xPos = coords[particle].transpose()[0]
    yPos = coords[particle].transpose()[1]
    zPos = coords[particle].transpose()[2]
    ax.scatter(xPos, yPos, zPos, zdir='z', s=20, c='k', depthshade=True)
    ax.set_xlim(0, 1024)
    ax.set_ylim(0, 1024)
    ax.set_zlim(0, 1024)

    ax.set_xlabel('X ')
    ax.set_ylabel('Y ')
    ax.set_zlabel('Z ')
    plt.show()
'''
