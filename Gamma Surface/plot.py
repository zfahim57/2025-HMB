import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import seaborn as sns
sns.set(font_scale=1.2, font='Arial')

def lammps_cell(xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz):
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo
    a = lx
    b = np.sqrt(ly**2 + xy**2)
    c = np.sqrt(lz**2 + xz**2 + yz**2)
    b = 10e-6 if b==0 else b
    c = 10e-6 if c==0 else c
    alpha = np.arccos((xy*xz + ly*yz)/(b*c))*180/np.pi
    beta = np.arccos(xz/c)*180/np.pi
    gamma = np.arccos(xy/b)*180/np.pi
#    print(a, b, c, alpha, beta, gamma)
    return a, b, c, alpha, beta, gamma

xyz = np.loadtxt('dump.minimized', skiprows=5, max_rows=3)
#print(xyz.shape)
xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz = xyz[0][1], xyz[0][0], xyz[1][1], xyz[1][0],xyz[2][1], xyz[2][0], xyz[0][2], xyz[2][2], xyz[1][2]
xlo -= min([0, xy, xz, xy+xz])
xhi -= max([0, xy, xz, xy+xz])
ylo -= min([0, yz])
yhi -= max([0, yz])

cell_par = lammps_cell(xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz) 
A = cell_par[0]*cell_par[1]*np.sin(cell_par[5])
#lx = xyz[0,1]-xyz[0,0]
#ly = xyz[1,1]-xyz[1,0]
#A = lx*ly
print(A)

xyz = np.loadtxt('logEf.txt')
xyz = xyz.T
#print(xyz.shape)
x = xyz[0]
y = xyz[1]
energy = xyz[3]

lcount = 101
#print(lcount)
xi = np.linspace(-5.3, 5.3, lcount)
yi = np.linspace(-6.3, 6.3, lcount)
xi, yi = np.meshgrid(xi, yi)

energy -= energy[len(energy)//2]
energy /= A
#print( energy[0])
energy = energy.reshape(lcount, lcount)
energy = energy.T
#ax_size = 13
#label_size = 16 
#cbar_size = 16
#print(energy.shape)
plt.contourf(xi, yi, energy, levels = 1000, cmap = 'viridis')
#cbar = plt.colorbar(label=r'Energy(kCal/mol-$\AA^2$)', format=tkr.FormatStrFormatter('%.3f'))
cbar = plt.colorbar(format=tkr.FormatStrFormatter('%.3f'))
cbar.set_label(r'Energy(kCal/mol-$\AA^2$)')#, fontsize=label_size) 
cbar.ax.tick_params()#labelsize=ax_size)
contours = plt.contour(xi, yi, energy, levels = 8, colors = 'black', linewidths=0.5)
plt.xlim(-3,3)
plt.ylim(-4,4)
plt.xlabel(r'X($\AA$)')#, fontsize=label_size)
plt.ylabel(r'Y($\AA$)')#, fontsize=label_size)
plt.tick_params(axis='x')#, labelsize=ax_size)
plt.tick_params(axis='y')#, labelsize=ax_size)
plt.tight_layout()
plt.savefig('Countour-00.jpg', dpi=500)


