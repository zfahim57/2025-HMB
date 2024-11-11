import matplotlib.pyplot as plt
import elastic
import mpl_toolkits.mplot3d.axes3d as axes3d
import numpy as np
import pandas as pd
import math
import matplotlib.colors as colors
import matplotlib as mpl
from matplotlib.ticker import LinearLocator
from matplotlib import cm
import seaborn as sns
import os
sns.set(font_scale=1.4, font='arial', style='white')

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


def spherical_grid(npoints=100):
    theta = np.linspace(0, np.pi, npoints)
    phi = np.linspace(0, 2 * np.pi, npoints)
    return np.meshgrid(theta, phi)
    # return theta, phi


def spherical_coord(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

def minimum_shear(mat, npoints=200):
    shear_val, t1, p1, c1 = 100, 0, 0, 0
    angle = np.linspace(0, 2*np.pi, npoints)

    for theta in angle:
        for phi in angle:
            #phi = 0
            for chi in angle:
                val = mat.shear([theta, phi, chi])
                if shear_val > val:
                    shear_val = val
                    t1, p1, c1 = theta, phi, chi
    return shear_val, t1, p1, c1

C = np.zeros((12, 6, 6))
T = 300
N = 50
shear_modulus = np.zeros((12, N))
chi = np.linspace(0, 2*np.pi, N)

data = np.loadtxt('graph/%g.txt' % T)
material = elastic.Elastic(data)
f = np.vectorize(lambda x, y: material.shear2D([x, y]))
theta, phi = spherical_grid(N)
r = f(theta, phi)
xmax, ymax, zmax = spherical_coord(r[1], theta, phi)
xmin, ymin, zmin = spherical_coord(r[0], theta, phi)
min_shear = np.min(r[0])
value, tht1, ph1, ch1 = minimum_shear(material, npoints=N)
tht1, ph1, ch1 = tht1*180/np.pi,360- ph1*180/np.pi, ch1*180/np.pi
#print(value, tht1, ph1, ch1)
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot( projection='3d')

norm = plt.Normalize(vmin=np.min(r[0]), vmax=np.max(r[0]))  # Normalize r between its min and max values
colors = cm.viridis(norm(r[0])) 

surf = ax.plot_surface(xmin, ymin, zmin, facecolors=colors, rstride=1, cstride=1, linewidth=0)

m = cm.ScalarMappable(cmap='viridis', norm=norm)
m.set_array(r[0])  # Set the array to be used for the colorbar
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
# Show the plot
#ax.text(2, -0.05, -0.15, f'$\\G_\\text{min}$ = {min_shear:0.02f} GPa \nat ($\\theta$, $\\phi$, $\\chi$) =  \n (${tht1:0.0f}$, ${ph1:0.0f}$, ${ch1:0.0f}$)', color='black', ha='center')
#ax.text(2, -0.05, -0.15, f'$G_{{min}}$ = {min_shear:0.02f} GPa at ($\\theta$, $\\phi$, $\\chi$) \n (${tht1:0.02f}$, ${ph1:0.02f}$, ${ch1:0.02f}$)', color='black', ha='center')
ax.text(
    2, -0.05, -0.15,
    f'$G_{{\\text{{min}}}} (\\theta, \\phi, \\chi)$\n'
    f'= {min_shear:.2f} GPa\n'
    f'at ({tht1:.0f}, {ph1:.0f}, {ch1:.0f})',
    color='black', ha='center'
)

ax.view_init(10,-60)
ax.set_axis_off()

o = -1.05
origin = [o, o, o]

length = 0.000001

ax.plot([origin[0], length], [origin[1], origin[1]], [origin[2], origin[2]], color='r', linewidth=2, label='X-axis')
ax.plot([origin[0], origin[0]], [origin[1], length], [origin[2], origin[2]], color='g', linewidth=2, label='Y-axis')
ax.plot([origin[0], origin[0]], [origin[1], origin[1]], [origin[2], length], color='b', linewidth=2, label='Z-axis')

ax.text(length+0.01, o, o, "X", color='red')#, fontsize=14)
ax.text(o, length+0.01, o, "Y", color='green')#, fontsize=14)
ax.text(o, o, length+0.01, "Z", color='blue')#, fontsize=14)

ax.dist = -100

cax = ax.inset_axes([0.27, 0.15, 0.5, 0.04])
cbar = fig.colorbar(m, ax=ax,orientation = 'horizontal', cax=cax)
cbar.set_label('Shear Modulus $G$ (GPa)')#, fontsize=18)
cbar.ax.tick_params()#labelsize=12)
move_ax(ax, dx=-0.2, axisoff=False)

plt.savefig('Fig6-minimum_shear.pdf')

