import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set(font_scale=1.6, font='Arial')

y_xz = [None] * 3

y_xy = np.load('300_xy.npy')*0.0433634*10e3
y_xz[0] = np.load('300_xz.npy')*0.0433634*10e3
y_yz = np.load('300_yz.npy')*0.0433634*10e3
y_xz[1] = np.load('250_xz.npy')*0.0433634*10e3
y_xz[2] = np.load('350_xz.npy')*0.0433634*10e3

#print(y_xy.shape)
x = []
srate = 2e-7
for i in range(1000):
    x.append(i*srate*10e2)

fig, axs = plt.subplots(2, 1, figsize=(10, 8))
axs[0].set_title('(a)')
axs[0].plot(x, y_xy[:1000], label='shear on xy plane', color='r')
axs[0].plot(x, y_xz[0][:1000], label='shear on xz plane', color='g')
axs[0].plot(x, y_yz[:1000], label='shear on yz plane', color='b')
axs[0].set_ylabel("Energy (meV)")
axs[0].legend(loc='upper right')
axs[0].grid(False)
axs[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

y_xz[0] -= y_xz[0][0]
y_xz[1] -= y_xz[1][0]
y_xz[2] -= y_xz[2][0]

axs[1].set_title('(b)')
axs[1].plot(x, y_xz[1][:1000], label='250K', color='darkgreen')
axs[1].plot(x, y_xz[0][:1000], label='300K', color='limegreen')
axs[1].plot(x, y_xz[2][:1000], label='330K', color='springgreen')
axs[1].set_ylabel("Energy (meV)")
axs[1].set_xlabel(r"Strain ($\epsilon$)")
axs[1].legend(loc='upper right')
axs[1].grid(False)

plt.tight_layout()
plt.savefig("Fig9-strain.pdf')
