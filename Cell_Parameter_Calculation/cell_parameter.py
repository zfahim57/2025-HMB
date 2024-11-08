from ase.io.lammpsrun import construct_cell
from ase.cell import Cell
import numpy as np

def cell_par(xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz): 
    xlo += min([0, xy, xz, xy+xz])
    xhi -= max([0, xy, xz, xy+xz])
    ylo += min([0, yz])
    yhi -= max([0, yz])
    diagdisp = [xlo, xhi, ylo, yhi, zlo, zhi]
    offdiag = [xy, xz, yz]
    cell, celldisp = construct_cell(diagdisp, offdiag)
    cell_cl = Cell.new(cell)
#    print(cell, celldisp)
    print(cell_cl.cellpar())
    return cell_cl.cellpar()

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
    print(a, b, c, alpha, beta, gamma)    
    return a, b, c, alpha, beta, gamma

#Initial Structure
xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz = 185.0826, 0, 88.112, 0, 172.5444, 0, 0.1165, -18.081, -2.9061 # This is from log file from lammps just after the initialization
cell_cl = cell_par(xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz)
lammps_cell(xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz)

#after MD equilibration
data = np.loadtxt('293_equilibrium.1000000.dump', skiprows=5, max_rows=3)

xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz = data[0][1], data[0][0], data[1][1], data[1][0],data[2][1], data[2][0], data[0][2], data[2][2], data[1][2]
#print('xhi',xhi,'xlo', xlo,'yhi', yhi,'ylo', ylo,'zhi', zhi,'zlo', zlo,'xy', xy,'yz', yz,'xz', xz)

xlo -= min([0, xy, xz, xy+xz])
xhi -= max([0, xy, xz, xy+xz])
ylo -= min([0, yz])
yhi -= max([0, yz])
cell_cl = cell_par(xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz)
cell_cl = lammps_cell(xhi, xlo, yhi, ylo, zhi, zlo, xy, yz, xz)

new_cell = Cell.fromcellpar([cell_cl[0],cell_cl[1],cell_cl[2],cell_cl[3],cell_cl[4],cell_cl[5]])
matrix = np.array([[1, 1, 2], [1, -1, 0], [4, 3, -1]])
inverse_matrix = np.linalg.inv(matrix)
im = np.linalg.inv([[6,0,0],[0,5,0],[0,0,3]])
temp = 0.5*np.eye(3)
original = inverse_matrix@(im@temp@new_cell)
original_cell = Cell.new(original)
print(original_cell.cellpar())
