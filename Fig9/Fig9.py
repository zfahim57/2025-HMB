import lammps_logfile
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import make_interp_spline
from scipy.signal import savgol_filter
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
sns.set(font_scale=1.6, font='arial')

from matplotlib.transforms import Transform
from matplotlib.ticker import FuncFormatter, FixedLocator
from matplotlib.scale import ScaleBase
from matplotlib.scale import register_scale
from matplotlib import scale


class PowerTransform(Transform):
    """
    A transformation that applies the square root scale.
    """
    input_dims = 1
    output_dims = 1
    is_separable = True

    def __init__(self, power=2):
        super().__init__()
        self.power = power

    def transform_non_affine(self, values):
        values = np.abs(values)
        result = np.float_power(values, 1/self.power)
        if np.any(np.isnan(result)):
            print(values)
            import sys; sys.exit()
        return np.float_power(values, 1/self.power)

    def inverted(self):
        return PowerInverseTransform(self.power)


class PowerInverseTransform(Transform):
    """
    Inverse of the square root transform.
    """
    input_dims = 1
    output_dims = 1
    is_separable = True

    def __init__(self, power=2):
        super().__init__()
        self.power = power

    def transform_non_affine(self, values):
        values = np.abs(values)
        result = np.float_power(values, self.power)
        if np.any(np.isnan(result)):
            print(values)
            import sys; sys.exit()
        return np.float_power(values, self.power)

    def inverted(self):
        return PowerTransform(self.power)

class PowerScale(ScaleBase):
    """
    A custom scale for a square root transformation.
    """
    name = 'power'

    def __init__(self, axis, power=2, **kwargs):
        super().__init__(axis)
        self.power = power

    def get_transform(self):
        return PowerTransform(self.power)

    def set_default_locators_and_formatters(self, axis):
        base = np.array([0, 100, 400, 1000, 2000, 4000])
        axis.set_major_locator(FixedLocator(base)) #np.float_power(base, self.power)))
        axis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.0f}"))

# Register the custom scale with Matplotlib
scale.register_scale(PowerScale)

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
p_300 = {}
pxz = {}
energy = {}
for var in direction:
    if var != '-450':
        log = lammps_logfile.File('%g_%s.log' % (300, var))
        d = 'P'+var
        z = abs(log.get(d, run_num=2))
        p_300[var] = z
        energy[var] = log.get('TotEng', run_num=2)*0.0433634*1e3
    else:
        log = lammps_logfile.File('300_-450.log')
        pxz = abs(log.get('Pxz', run_num = 2))
        pyz = abs(log.get('Pyz', run_num = 2))
        z = abs(pyz*np.sin(theta) - pxz*np.cos(theta))
        p_300[var] = z
        energy[var] = log.get('TotEng', run_num=2)*0.0433634*1e3

fig, axs = plt.subplots(2, 1, figsize=(9, 8))
for label in ['xy', 'xz', 'yz', '-450']:
    y = p_300[label]
    y = savgol_filter(y, 201, 2)
    y -= y[0] + 1e-2
    e = savgol_filter(energy[label], 201, 2)
    if label == 'xy':
        id = int(len(y)/2.5)
    else:
        id = len(y) - 1

    if label == '-450':
        x0 = xp.copy()
        #label = f'$\sigma_{{[\\overline{{4}}50](001)}}$'
        label = f'$\sigma_{{\\text{{min}}}}$'
    else:
        label = '$\sigma_{' + label + '}$'
        x0 = x.copy()
    axs[0].plot(x0[:id], y[:id], lw=2, label = label)
    axs[1].plot(x0[:id], e[:id], lw=2, label = label)

#axs[0].set_yscale('power', power=3)
axs[0].set_yscale('power', power=2.3)
axs[0].set_title('(a)')
axs[1].set_title('(b)')
axs[0].set_xticklabels([])
#axs[0].grid(False)
axs[1].grid(False)
axs[0].set_xlim(-0.001, 0.2)
axs[1].set_xlim(-0.001, 0.2)
axs[0].set_ylim(0, 4200)
#axs[1].set_ylim(0, 690)
axs[0].legend(loc=1)
axs[1].legend(loc=1)
axs[0].set_ylabel('$\sigma$ (MPa)')
axs[1].set_ylabel('$E$ (meV)')
axs[1].set_xlabel('$\epsilon$')
axs[1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
plt.tight_layout()
plt.savefig('Fig9-strain.pdf')
plt.show()
