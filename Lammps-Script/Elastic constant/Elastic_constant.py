import numpy             as np
import scipy.stats       as stats
import matplotlib        as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker   import FormatStrFormatter


# Variables

Temp  = [300]   # Temperatures
T     = len(Temp)
tnpt  = [500]   # NPT equilibration time
tequi = 8000               # Equilibration time
delta = 0.001              # Small deformation
C     = np.zeros((T,6,6))  # Elastic Constants
C_err = np.zeros((T,6,6))  # Std dev. of Elast. Const.

C11_all = np.zeros(T)      # C11 as a mean of c11, c22 and c33
C12_all = np.zeros(T)      # C12 as a mean of c12, c21, c13, c31, c23 and c32
C44_all = np.zeros(T)      # C44 as a mean of c44, c55 and c66
C11_allerr = np.zeros(T)      # C11_all err
C12_allerr = np.zeros(T)      # C12 all err
C44_allerr = np.zeros(T)      # C44 all err
B       = np.zeros(T)      # bulk modulus
G       = np.zeros(T)      # shear modulus
mu      = np.zeros(T)      # Poisson’s ratio

# method 2

C_v2     = np.zeros((T,6,6))  # Elastic Constants
C_err_v2 = np.zeros((T,6,6))  # Std dev. of Elast. Const.
C11_all_v2 = np.zeros(T)      # C11 as a mean of c11, c22 and c33
C12_all_v2 = np.zeros(T)      # C12 as a mean of c12, c21, c13, c31, c23 and c32
C44_all_v2 = np.zeros(T)      # C44 as a mean of c44, c55 and c66
C11_allerr_v2 = np.zeros(T)      # C11_all err
C12_allerr_v2 = np.zeros(T)      # C12 all err
C44_allerr_v2 = np.zeros(T)      # C44 all err
B_v2       = np.zeros(T)      # bulk modulus
G_v2       = np.zeros(T)      # shear modulus
mu_v2      = np.zeros(T)      # Poisson’s ratio

l1 = np.zeros(6)
l2 = np.zeros(6)
s1 = np.zeros(6)
s2 = np.zeros(6)

# Operations

for i in np.arange(T):
    
    for j in np.arange(6):
        data2 = np.loadtxt('equi_nvtneg_%gK_%g.dat' % (Temp[i],j+1), skiprows=50)
        data3 = np.loadtxt('equi_nvtpos_%gK_%g.dat' % (Temp[i],j+1), skiprows=50)
                
        l1[0] = np.array(data2[:,3]).mean()
        l1[1] = np.array(data2[:,4]).mean()
        l1[2] = np.array(data2[:,5]).mean()
        l1[3] = np.array(data2[:,8]).mean()
        l1[4] = np.array(data2[:,7]).mean()
        l1[5] = np.array(data2[:,6]).mean()
        
        l2[0] = np.array(data3[:,3]).mean()
        l2[1] = np.array(data3[:,4]).mean()
        l2[2] = np.array(data3[:,5]).mean()
        l2[3] = np.array(data3[:,8]).mean()
        l2[4] = np.array(data3[:,7]).mean()
        l2[5] = np.array(data3[:,6]).mean()
        
        s1[0] = np.array(data2[:,3]).std()
        s1[1] = np.array(data2[:,4]).std()
        s1[2] = np.array(data2[:,5]).std()
        s1[3] = np.array(data2[:,8]).std()
        s1[4] = np.array(data2[:,7]).std()
        s1[5] = np.array(data2[:,6]).std()
        
        s2[0] = np.array(data3[:,3]).std()
        s2[1] = np.array(data3[:,4]).std()
        s2[2] = np.array(data3[:,5]).std()
        s2[3] = np.array(data3[:,8]).std()
        s2[4] = np.array(data3[:,7]).std()
        s2[5] = np.array(data3[:,6]).std()
        
        for k in np.arange(6):
            C_v2[i][j][k] = (l1[k]-l2[k])/0.02*(1.01325e-4)
            C_err_v2[i][j][k] = np.sqrt(s1[k]**2+s2[k]**2)*(1.01325e-4)
            
    
    # method 2
    
    C11_all_v2[i] = np.array([C_v2[i][0][0],C_v2[i][1][1],C_v2[i][2][2]]).mean()
    C12_all_v2[i] = np.array([C_v2[i][0][1],C_v2[i][0][2],C_v2[i][1][0],C_v2[i][1][2],C_v2[i][2][0],C_v2[i][2][1]]).mean()
    C44_all_v2[i] = np.array([C_v2[i][3][3],C_v2[i][4][4],C_v2[i][5][5]]).mean()
    
    C11_allerr_v2[i] = np.sqrt((np.array([C_v2[i][0][0],C_v2[i][1][1],C_v2[i][2][2]]).std())**2/3+max(C_err_v2[i][0][0],C_err_v2[i][1][1],C_err_v2[i][2][2])**2)
    C12_allerr_v2[i] = np.sqrt((np.array([C_v2[i][0][1],C_v2[i][0][2],C_v2[i][1][0],C_v2[i][1][2],C_v2[i][2][0],C_v2[i][2][1]]).std())**2/6+max(C_err_v2[i][0][1],C_err_v2[i][0][2],C_err_v2[i][1][0],C_err_v2[i][1][2],C_err_v2[i][2][0],C_err_v2[i][2][1])**2)
    C44_allerr_v2[i] = np.sqrt((np.array([C_v2[i][3][3],C_v2[i][4][4],C_v2[i][5][5]]).std())**2/3+max(C_err_v2[i][3][3],C_err_v2[i][4][4],C_err_v2[i][5][5])**2)

    print("Temperature ", Temp[i])
    
    print("method 2")
    np.set_printoptions(precision=3,suppress=True)
    print("C_v2 = ", C_v2[i])
    
    print("C_err_v2 = ", C_err_v2[i])
    
    np.set_printoptions(precision=3,suppress=True)
    print("C11_all_v2 = ", C11_all_v2[i], " +- ", C11_allerr_v2[i], " GPa")
    print("C12_all_v2 = ", C12_all_v2[i], " +- ", C12_allerr_v2[i], " GPa")
    print("C44_all_v2 = ", C44_all_v2[i], " +- ", C44_allerr_v2[i], " GPa")

    Cn = np.zeros_like(C_v2)
    Cn = C_v2 + C_err_v2

    print(Cn)

new_file = open("Elastic_constant.txt", "w")
cw = 25
for i in range(6):
    for j in range(6):
        if i == j:
            new_file.write(f'{C_v2[0,i,i]:<{cw}}')

        else:
            avg = (C_v2[0,i,j]+C_v2[0,j,i])/2
            new_file.write(f'{avg:<{cw}}')
    new_file.write('\n')

