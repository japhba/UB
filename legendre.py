from numpy.polynomial import Legendre
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns
import itertools

palette1 = itertools.cycle(sns.color_palette("Blues_r"))
palette2 = itertools.cycle(sns.color_palette("Blues_r"))

x = np.linspace(-1, 1, 1000)
y = np.linspace(-np.pi, np.pi, 1000)

f0 = Legendre((1, 0, 0))(np.cos(y))
f1 = Legendre((0, 1, 0))(np.cos(y))
f2 = Legendre((0, 0, 1))(np.cos(y))
f3 = Legendre((0, 0, 0, 1))(np.cos(y))
f4 = Legendre((0, 0, 0, 0, 1))(np.cos(y))

f = [f0, f1, f2, f3, f4]
F = []


for k in range(5):
    FF = f0*0
    for i in range(k+1):
        for j in range(k+1):
            FF += f[i]*f[j]
    F.append(FF)

# figure setup
fig = plt.figure(figsize = (10,3))
gs  = fig.add_gridspec(1, 3)
ax1 = fig.add_subplot(gs[0, 0:2])
ax2 = fig.add_subplot(gs[0, 2], projection='polar')
xT  = ax2.get_xticks()
xL  = ['0',r'$45^\circ$',r'$90^\circ$',r'$135^\circ$',
    r'$\pm180^\circ$',r'$-135^\circ$', r'$-90^\circ$', r'$-45^\circ$']

ax2.set_xticks(xT)
ax2.set_xticklabels(xL)

ax2.set_yticklabels([])

[ax1.plot(y, np.abs(f)**2, c = next(palette1), label = "$P_" + str(i) + "$") for f,i in zip(f, range(len(f)))]
[ax2.plot(y, np.abs(f)**2, c = next(palette2)) for f,i in zip(f, range(len(f)))]

xT  = ax1.get_xticks()
w = np.array([-180,-135,-90,-45,0,45,90,135,180])
xL = [str(xTT) + "$^\\circ$" for xTT in w]

ax1.set_xticks(w*2*np.pi/360)
ax1.set_xticklabels(xL)

#plt.subplots_adjust(right = 1, left = 1)

ax1.legend()

plt.savefig("legendre.pdf")
