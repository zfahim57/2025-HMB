import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
from matplotlib.ticker import ScalarFormatter
import seaborn as sns
sns.set(font_scale=1.5, font='Arial')

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
print(A)

xyz = np.loadtxt('logEf.txt')
xyz = xyz.T
x = xyz[0]
y = xyz[1]
energy = xyz[3]

lcount = 101
xi = np.linspace(-5.3, 5.3, lcount)
yi = np.linspace(-6.3, 6.3, lcount)
xi, yi = np.meshgrid(xi, yi)

energy -= energy[len(energy)//2]
energy /= (A*15360)
energy = energy.reshape(lcount, lcount)
energy = energy.T
vmin = 0.0
vmax = 0.041
energy[energy > vmax] = vmax
energy[energy < vmin] = vmin
energy *= 0.043363 * 1000; print(14.38 * 0.043363)

contour = plt.contourf(xi, yi, energy, levels=25, cmap='viridis')
#cbar = plt.colorbar(contour, format=tkr.FormatStrFormatter('%.6f'))
cbar = plt.colorbar(contour, format=tkr.FuncFormatter(lambda x, pos: f"{x * 1e4:.2f}"))
cbar.set_label(r'$\Delta E$ (meV / atom-$\mathrm{\AA}^2$)')
cbar.ax.tick_params()
cbar.ax.text(
    1.8, 1.02,  # Position: slightly above the top of the bar
    r"($\times 10^{-4}$)",  # LaTeX formatted scaling factor
    transform=cbar.ax.transAxes,  # Use colorbar's coordinate system
    ha='center',  # Horizontal alignment
    va='bottom',  # Vertical alignment
    fontsize=10
)
# Coordinates of points O, A, and B
x_O, y_O = -0.2, 0.05
x_A, y_A = -1.6, 1.15
x_B, y_B = 1.9, 0.6

# Plot points O, A, and B
plt.text(x_O-0.1, y_O-0.3, "O", color='white')
plt.text(x_A-0.2, y_A-0.3, "A", color='white')
plt.text(x_B+0.05, y_B-0.3, "B", color='white')

# Draw arrows from O to A and O to B
plt.arrow(x_O, y_O, x_A - x_O, y_A - y_O, head_width=0.1, head_length=0.1, fc='white', ec='white')
plt.arrow(x_O, y_O, x_B - x_O, y_B - y_O, head_width=0.1, head_length=0.1, fc='white', ec='white')

plt.xlim(-2.8, 2.8)
plt.ylim(-2.8, 2.8)
plt.xlabel(r'$X$ ($\mathrm{\AA}$)')
plt.ylabel(r'$Y$ ($\mathrm{\AA}$)')
plt.tick_params(axis='x')#, labelsize=ax_size)
plt.tick_params(axis='y')#, labelsize=ax_size)
plt.tight_layout()
plt.savefig('Fig5-gamma_surface.pdf')

