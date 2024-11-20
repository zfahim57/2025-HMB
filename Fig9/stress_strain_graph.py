import lammps_logfile
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import make_interp_spline
from scipy.signal import savgol_filter
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
sns.set(font_scale=1.8, font='arial')

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
xyz = np.loadtxt('equilibrium.40000.dump', skiprows=5, max_rows=3)
#print(xyz.shape)
xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz = xyz[0][1], xyz[0][0], xyz[1][1], xyz[1][0],xyz[2][1], xyz[2][0], xyz[0][2], xyz[2][2], xyz[1][2]
xlo -= min([0, xy, xz, xy+xz])
xhi -= max([0, xy, xz, xy+xz])
ylo -= min([0, yz])
yhi -= max([0, yz])

cell_par = lammps_cell(xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz)
lx, ly, lz = cell_par[0], cell_par[1], cell_par[2]
a, b = 1, 1.25
c = np.sqrt(a**2 + b**2)
theta = np.arctan((b*ly)/(a*lx))

log1 = lammps_logfile.File("300_xz.log")

lz1 = log1.get("Lz", run_num=2)
x = []
srate = 2e-7
#print(len(lz1))
for i in range(len(lz1)):
    x.append(i*srate*10e2)
xp = []
for i in range(len(lz1)):
    xp.append(i*srate*10e2*c)

direction = ['xy', 'xz','yz', '-450']
#d_var = ['xy', 'xz', 'yz', '-110']
filter_param = 201
p_300 = {}
pxz = {}
for var in direction:
    if var != '-450':
        log = lammps_logfile.File('%g_%s.log' % (300, var))
        d = 'P'+var
        z = abs(log.get(d, run_num=2))
        p_300[var] = z
    else:
        log = lammps_logfile.File('300_-450.log')
        pxz = abs(log.get('Pxz', run_num = 2))
        pyz = abs(log.get('Pyz', run_num = 2))
        z = abs(pyz*np.sin(theta) - pxz*np.cos(theta))
        p_300['-450'] = z
        print(p_300['-450'])
#print(np.mean(p_300['-110']))
fig, axs = plt.subplots(2, 1, figsize=(9, 8))

for i in range(4):
    var = direction[i]
    if(i != 3):
        #d = var[1:3]
        y = p_300[var]
        yf = savgol_filter(y, filter_param, 2)
        axs[0].plot(x, yf, label = f'$\sigma_{{{var}}}$')
#            axs[0].legend(True)
    else:
        y = p_300[var]
        #print(y.shape)
        y1 = savgol_filter(y, filter_param, 2)
        y = p_300['xz']
        y2 = savgol_filter(y, filter_param, 2)
        axs[1].plot(xp, y1, label = f'$\sigma_{{[\\overline{{4}}50](001)}}$')
        axs[1].plot(x, y2, label = f'$\sigma_{{xz}}$')
     
axs[0].set_title('(a)')
axs[1].set_title('(b)')
axs[0].set_xticklabels([])
axs[0].tick_params('y')
axs[1].tick_params('y')
axs[1].tick_params('x')
axs[0].grid(False)
axs[1].grid(False)
[x.set_xlim(0, 0.2) for x in axs]
axs[0].legend(loc = 'center right', fontsize=14)
axs[1].legend(loc = 'center right', fontsize=14)
axs[0].set_ylabel('$\sigma$ (MPa)')
axs[1].set_ylabel('$\sigma$ (MPa)')
axs[1].set_xlabel('$\epsilon$')
axs[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.savefig('Fig9-strain.pdf')
plt.show()
