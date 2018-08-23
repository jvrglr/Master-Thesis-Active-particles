import numpy as np
import matplotlib.pyplot as plt
import time
import pylab

R=2.0
for i in range(1,501):
    force=np.loadtxt("force."+str(i)+".dat")
    #self=np.loadtxt("self_U."+str(i)+".dat")
    radial= np.loadtxt("radial."+str(i)+".dat")
    my_dpi=96
    #fig = plt.figure(figsize=(1600/my_dpi, 1200/my_dpi), dpi=my_dpi)
    fig = plt.figure()
    #ax = fig.add_axes([-0.38, -0.38, 0.9, 0.9]) # main axes
    #ax2 = fig.add_axes([0.1, 0.1, 0.4, 0.4]) # main axes
    ax = fig.add_axes([0.1, 0.1, 0.9, 0.9]) # main axes
    ax2 = fig.add_axes([0.57, 0.57, 0.4, 0.4]) # main axes
    #ax.axis('off')
    q=ax.quiver(force[:,0], force[:,1], force[:,2], force[:,3], color="midnightblue",angles='xy', scale_units='x',scale=20,width=0.001)
    m=ax.quiver(force[:,0], force[:,1], force[:,4], force[:,5], color="red",angles='xy', scale_units='x',scale=20,width=0.001)
    radio=ax.quiver(0.2, 0.2, R, 0.0, color="black",angles='xy', scale_units='x',scale=1,width=0.005)

    ax2.plot(radial[:,0],radial[:,1],"-")
    ax2.quiverkey(m,X=0.75,Y=0.90,U=8,label='Self propulsion', labelpos='E')
    ax2.quiverkey(q,X=0.75,Y=0.85,U=8,label='Repulsive force', labelpos='E')

    fig.savefig("image."+str(i)+".png", bbox_inches='tight',dpi=300)
    plt.close(fig)
