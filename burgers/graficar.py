import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot

x = np.linspace(0,2,41)
y = np.linspace(0,2,41)

for i in range (1,501):
    u1 = np.loadtxt("matrizdata"+str(i)+".dat")
    d1 = np.split(u1,41)
    u = np.asarray(d1)
    
    fig = pyplot.figure(figsize=(11,7))
    ax = pyplot.gca(projection='3d')
    X,Y = np.meshgrid(x,y)
    wire1 = ax.plot_wireframe(X,Y,u[:], cmap=cm.coolwarm)
    #wire1 = ax.plot_wireframe(X,Y,v[:], cmap=cm.coolwarm)
    
    plt.savefig('burger'+str(i)+'.png')
    plt.close(fig)
os.system("convert -delay 1 $(ls -v burger*.png) burger.gif")
