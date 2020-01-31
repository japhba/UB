import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ame import *
from plotinit import initMPL
#initMPL()

Z, N = np.mgrid[1:118+1, 1:177+1]
A = Z + N
E_bethe = bethe(Z,N)


try:
    E_ame = np.load("ameGrid.npy")
except:
    E_ame = binding(Z,N, nominal = True)
    np.save("ameGrid.npy", E_ame)

E_ame /= 1000 # ame data in MeV
E_ame[E_ame < 1e-2] = np.nan

fig = plt.figure(figsize = (12,9))
ax1, ax2, ax3 , ax4 = plt.subplot(221), plt.subplot(222), plt.subplot(223), plt.subplot(224)

im = ax1.matshow(E_bethe/A, origin = "lower", extent=[1-0.5,178-0.5,1-0.5,119-0.5])
ax1.set_title("Semi empirical mass formula")
ax1.xaxis.set_ticks_position('bottom')
ax1.set_ylabel("proton number $Z$")
ax1.set_xlabel("neutron number $N$")

divider = make_axes_locatable(ax1)
cax1 = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax1, orientation='vertical')


im = ax2.matshow(E_ame, origin = "lower", extent=[1-0.5,178-0.5,1-0.5,119-0.5])
ax2.set_title("Experimental values (AME2016)")
ax2.xaxis.set_ticks_position('bottom')
ax2.set_ylabel("proton number $Z$")
ax2.set_xlabel("neutron number $N$")

divider = make_axes_locatable(ax2)
cax2 = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im, cax=cax2, orientation='vertical')

im = ax3.matshow(-(E_bethe/A - E_ame)*1000, origin = "lower", cmap = "viridis", vmin = -130, vmax = 130, extent=[1-0.5,178-0.5,1-0.5,119-0.5])
ax3.set_title("Difference")
ax3.xaxis.set_ticks_position('bottom')
ax3.set_ylabel("proton number $Z$")
ax3.set_xlabel("neutron number $N$")

divider = make_axes_locatable(ax3)
cax3 = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(im, cax=cax3, orientation='vertical')
cb.ax.set_ylabel("$E_B$ in keV")

l = 5
E_ame2 = E_ame[92 - l-1:92 + l-1 , 238-92-l-1:238-92+l-1]
im = ax4.matshow((E_ame2-E_ame[92-1, 238-92-1])*1000, origin = "lower", extent=[238-92-l-.5,238-92+l-0.5 , 92 - l-0.5,92 + l-0.5], cmap = "bwr")
ax4.set_title("Relative energy to $^{238}$U")
ax4.xaxis.set_ticks_position('bottom')
ax4.set_ylabel("proton number $Z$")
ax4.set_xlabel("neutron number $N$")

# ax4.set_xlim(238-92-l, 238-92+l)
# ax4.set_ylim(92 - l, 92 + l)

ax4.scatter([238-92],[92], marker = "s", facecolors='none', edgecolors = "black", s = 70)

divider = make_axes_locatable(ax4)
cax2 = divider.append_axes('right', size='5%', pad=0.05)
cb = fig.colorbar(im, cax=cax2, orientation='vertical')
cb.ax.set_ylabel("$\\Delta E_B$ in keV")

# for ax in [ax1, ax2, ax3]:
#     labels = [item.get_text() for item in ax.get_xticklabels(which = 'both')]
#     print(labels)
#     ax.set_xticklabels(labels)

plt.tight_layout()
plt.savefig("bindings.pdf")

plt.show()
