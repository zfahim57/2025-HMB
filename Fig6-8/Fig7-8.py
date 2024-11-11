import elastic
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import pandas as pd
import math
from matplotlib.colors import LightSource
from matplotlib import cbook, cm
from matplotlib.axes import Axes
import seaborn as sns
sns.set(font_scale=1.5, font='Arial')

def move_ax(ax, dx=None, dy=None, axisoff=True):
    pos1 = ax.get_position()
    pts1 = pos1.get_points()

    if dx is not None:
        pts1[1][0] += dx
        pts1[0][0] += dx
    if dy is not None:
        pts1[1][1] += dy
        pts1[0][1] += dy

    pos1.set_points(pts1)
    ax.set_position(pos1)
    if axisoff:
        ax.axis('off')


def spherical_coord(r, theta):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y

def shear_calculation(Tstart, delT):
    T, dif = Tstart, delT
    n, m = (400-T)//dif + 1, 1000
    C = np.zeros((n, 6, 6))
    shear_modulus = np.zeros((n, m))
    chi = np.linspace(0, 2*np.pi, m)
    x = np.zeros((3, n, m))
    y = np.zeros((3, n, m))
    z = np.zeros((3, n, m))

    for i in range(n):
        theta = 0
        phi = 0
        data = np.loadtxt('graph/%g.txt' % T)
        mat = elastic.Elastic(data)
        # smin, smax, chi1, chi2 = mat.shear3D(0, 0)
        # print(smin, chi1*180/np.pi)
    #    print(mat.shear([theta, phi, 0]))
        r = []
        for c in chi:
            r.append(mat.shear([theta, phi, c]))
        shear_modulus[i] = r
        x[2, i], y[2, i] = spherical_coord(r, chi)
        z[2, i] = T  # np.sqrt(x[i]**2 + y[i]**2)

        theta = np.pi/2
        r = []
        for c in chi:
            r.append(mat.shear([theta, phi, c]))
        shear_modulus[i] = r
        x[0, i], y[0, i] = spherical_coord(r, chi)
        z[0, i] = T  # np.sqrt(x[i]**2 + y[i]**2)

        phi = np.pi/2
        r = []
        for c in chi:
            r.append(mat.shear([theta, phi, c]))
        shear_modulus[i] = r
        x[1, i], y[1, i] = spherical_coord(r, chi)
        z[1, i] = T  # np.sqrt(x[i]**2 + y[i]**2)

        T = T + dif
    return x, y, z

#plotting the xy plane from 350 to 400
x, y, z = shear_calculation(320, 10)

scatter = plt.scatter(x[2], y[2], c=z[2], s=5, cmap='viridis')
cbar = plt.colorbar(scatter)
cbar.set_label("Temperature (K)")#, fontsize=16)
cbar.set_ticks(np.linspace(320, 400, num=9))
#plt.plot(lx, ly, color='r', linestyle='->')
plt.xlabel('$G_X$ (GPa)')
plt.ylabel('$G_Y$ (GPa)')
ax = plt.gca()
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
# Coordinates of points O, A, and B
x_O, y_O = 0, 0
x_A, y_A = -0.26, 0.25

# Plot points O, A, and B
plt.text(x_A-0.18, y_A+0.1, r"[$\overline{1}10]$", color='k')

plt.arrow(x_O, y_O, x_A - x_O, y_A - y_O, head_width=0.05, head_length=0.1, fc='k', ec='k')
 
plt.tight_layout()
plt.grid(False)
plt.savefig('Fig8-hightemp_shear.pdf')
plt.close()

sns.set(font_scale=2.6, font='Arial')
#plotting xz and yz from 50 to 400
x, y, z = shear_calculation(50, 50)
#fig, axs = plt.subplots(1, 2, figsize=(12, 7.5))
fig, axs = plt.subplots(1, 2, figsize=(12, 7.5), gridspec_kw={'width_ratios': [1, 1.2]})
plt.subplots_adjust(wspace=0.25)
labelpad_value = -10
for i in range(2):
    scatter = axs[i].scatter(x[i], y[i], c=z[i], s=5, cmap='viridis')
    if i == 2:
        axs[i].set_xlabel('X')
        axs[i].set_ylabel('Y')
        axs[i].set_title('(c)')
        axs[i].tick_params(axis='both', which='major', pad=10)
    if i == 0:
       axs[i].set_xlabel('$G_x$ (GPa)')
       axs[i].set_ylabel('$G_z$ (GPa)')
       axs[i].set_title('(a)')
       axs[i].tick_params(axis='both', which='major', pad=10)
    if i == 1:
       axs[i].set_xlabel('$G_x$ (GPa)')
       #axs[i].set_ylabel('Z')
       axs[i].set_title('(b)')
       axs[i].tick_params(axis='both', which='major', pad=10)
    axs[i].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axs[i].grid(False)

move_ax(axs[0], dx=0.005, axisoff=False)
move_ax(axs[1], dx=-0.005, axisoff=False)

cbar = plt.colorbar(scatter)
cbar.set_label("Temperature (K)")#, fontsize=16)
cbar.set_ticks(np.linspace(50, 400, num=8))
fig.tight_layout()
plt.savefig('Fig7-shear_projection.pdf')
