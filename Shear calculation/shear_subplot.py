from ELATE import elastic
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import pandas as pd
import math
from matplotlib.colors import LightSource
from matplotlib import cbook, cm
from matplotlib.axes import Axes
import seaborn as sns
sns.set(font_scale=1.4, font='Arial')


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
x, y, z = shear_calculation(350,10)
lx = np.linspace(0, -0.2, 10)
angle = 2.1746504934615634
ly = np.tan(angle) * lx

scatter = plt.scatter(x[2], y[2], c=z[2], cmap='viridis')
cbar = plt.colorbar(scatter)
cbar.set_label("Temperature (K)")#, fontsize=16)
cbar.ax.tick_params()#labelsize=15)
plt.axhline(y=0, color='black', linestyle='-')
plt.axvline(x=0, color='black', linestyle='-')
plt.plot(lx, ly, color='r', linestyle='-')
plt.xlabel('X')#, fontdict={'family': 'serif', 'size': 16})
plt.ylabel('Y')#, fontdict={'family': 'serif', 'size': 16})
plt.tick_params(axis='x')#, labelsize=15)
plt.tick_params(axis='y')#, labelsize=15)
plt.tight_layout()
plt.grid(False)
plt.savefig('350_400_shear_line-1.jpg', dpi=500)
plt.show()
plt.close()

#plotting xz and yz from 50 to 400
x, y, z = shear_calculation(50, 50)
fig, axs = plt.subplots(1, 2, figsize=(12, 8))
plt.subplots_adjust(wspace=0.5)
label_size = 20
abc_size = 22
ax_size = 18
labelpad_value = -10
for i in range(2):
    scatter = axs[i].scatter(x[i], y[i], c=z[i], cmap='viridis')
    if i == 2:
        axs[i].set_xlabel('X')#, fontdict={'family': 'Arial', 'size': label_size})
        axs[i].set_ylabel('Y')#, fontdict={'family': 'Arial', 'size': label_size}, labelpad=labelpad_value)
        axs[i].set_title('(c)')#, fontdict={'family': 'Arial', 'size': abc_size})
        axs[i].tick_params(axis='both', which='major', pad=10)
    if i == 0:
       axs[i].set_xlabel('Y')#, fontdict={'family': 'Arial', 'size': label_size})
       axs[i].set_ylabel('Z')#, fontdict={'family': 'Arial', 'size': label_size}, labelpad=labelpad_value)
       axs[i].set_title('(a)')#, fontdict={'family': 'Arial', 'size': abc_size})
       axs[i].tick_params(axis='both', which='major', pad=10)
    if i == 1:
       axs[i].set_xlabel('X')#, fontdict={'family': 'Arial', 'size': label_size})
       axs[i].set_ylabel('Z')#, fontdict={'family': 'Arial', 'size': label_size}, labelpad=labelpad_value)
       axs[i].set_title('(b)')#, fontdict={'family': 'Arial', 'size': abc_size})
       axs[i].tick_params(axis='both', which='major', pad=10)
    axs[i].tick_params(axis='x')#, labelsize=ax_size)
    axs[i].tick_params(axis='y')#, labelsize=ax_size)
    axs[i].set_xlim(-3,3)
    axs[i].set_ylim(-12,12)
    axs[i].grid(False)
#        move_ax(axs[i], dx=-0.2, axisoff=False)
fig.subplots_adjust(left=0.35)
cbar = plt.colorbar(scatter)
cbar.set_label("Temperature (K)")#, fontsize=abc_size)
cbar.ax.tick_params(labelsize=ax_size) 
fig.tight_layout()
plt.subplots_adjust(hspace=0.6, wspace=0.5)
plt.savefig('subplot_2.jpg',dpi=500)
plt.show()



