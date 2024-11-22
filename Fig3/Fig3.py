import lammps_logfile
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import make_interp_spline
from scipy.signal import savgol_filter
from matplotlib.ticker import FuncFormatter
import seaborn as sns
sns.set(font_scale=2, font='arial')

log = lammps_logfile.File("log.lammps")

t = log.get("Temp", run_num = 2)
p_energy1 = log.get("PotEng", run_num = 2)*0.0433634*1000
lzz1 = log.get("Lz", run_num = 2)/48

t2 = log.get("Temp", run_num = 4)
p_energy2 = log.get("PotEng", run_num = 4)*0.0433634*1000
lzz2 = log.get("Lz", run_num = 4)/48
#p_energy2 = np.flip(p_energy2)
#lzz2 = np.flip(lzz2)

temp = savgol_filter(t,71,2)
temp2 = savgol_filter(t2,71,2)
Energy1 = savgol_filter(p_energy1,71,2)
lzf1 = savgol_filter(lzz1,71,2)

Energy2 = savgol_filter(p_energy2,71,2)
lzf2 = savgol_filter(lzz2,71,2)

linewidth = 2.5
formatter = FuncFormatter(lambda x, pos: f'{x:.2f}')
fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(10,8))
ax1.plot(temp,lzf1, label='Heating', color='brown', lw=linewidth)
ax1.plot(temp2,lzf2, label='Cooling', color='teal', lw=linewidth)
ax1.set_ylabel(r'$d$ ($\AA$)')
ax1.legend(loc='lower right')
ax1.set_title("(a)")

ax2.plot(temp,Energy1, label='Heating', color='brown', lw=linewidth)
ax2.plot(temp2,Energy2, label='Cooling', color='teal', lw=linewidth)
ax2.set_ylabel("$PE$ (meV / atoms)")
ax2.set_xlabel("Temperature(K)")
ax2.set_title("(b)")
ax1.tick_params(axis='y')
ax2.tick_params(axis='y')
ax2.tick_params(axis='x')
ax1.grid(False)
ax2.grid(False)
ax1.yaxis.set_major_formatter(formatter)
plt.savefig("Fig3-heatingcycle.pdf")
