import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

datos = sys.argv[1]
nombreSinDat = datos.strip('.dat')
arregloStrings = nombreSinDat.split('_')
x0 = float(arregloStrings[1])
y0 = float(arregloStrings[2].strip('.dat'))
data = np.loadtxt(datos)

fig2D = plt.figure()
fig2D.set_size_inches(18.5,10.5)
plt.xlabel('x(m)')
plt.ylabel('y(m)')
plt.title('Trayectoria en el plano xy con Energia '+str(x0)+' y pitch angle' + str(y0))
plt.plot(data[:,1],data[:,2])

fig3D = plt.figure()
fig3D.set_size_inches(18.5,10.5)
ax = fig3D.add_subplot(1, 2, 1, projection='3d')

u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 6378.1*1000 * np.outer(np.cos(u), np.sin(v))
y = 6378.1*1000 * np.outer(np.sin(u), np.sin(v))
z = 6378.1*1000 * np.outer(np.ones(np.size(u)), np.cos(v))

ax.set_xlabel('X(m)')
ax.set_ylabel('Y(m)')
ax.set_zlabel('Z(m)')
plt.title('Trayectoria en 3D con Energia' + str(x0)+' y pitch angle' + str(y0))
ax.plot_surface(x, y, z, color='b')
plt.plot(data[:,1],data[:,2],data[:,3], color='r' )

pp = PdfPages('t'+nombreSinDat+'.pdf')
pp.savefig(fig2D)
pp.savefig(fig3D)
pp.close()
