import numpy as np
import matplotlib.pyplot as plt

C = np.zeros((12,6,6))
T = 300
x = []
for i in range(12):
    data = np.loadtxt('%g.txt'% T)
    #print(data.dtype)
    print(data.shape)
    #print(data)
    C[i] = data
    x.append(T)
    T = T+10

#y1 = C[:,0,0]
y4 = C[:, 3,3]
y5 = C[:, 4,4]
#print(y)
plt.plot(x,y4, label='C44')
plt.plot(x, y5, label='C55')
plt.legend()
plt.savefig('C44_55.jpg', dpi=500)
