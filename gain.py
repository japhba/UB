from ame import excess, mass, bethe, binding
import matplotlib.pyplot as plt
import numpy as np

import pandas as pd

# plot energy gain graph
A = 236
Z = 92
N = A-Z


def yields(Z,N):
    df = pd.read_csv('/Users/jan/Google Drive/Uni/Uebungen/7. Semester/Nuclear Physics/yields2.csv', sep = ";")

    row = df.loc[(df['Z'] == Z) & (df['N'] == N)]
    try:
        return row["Yield"].to_numpy()[0]
    except:
        return np.NaN

yields = np.vectorize(yields, signature = '(),()->()')

def E2Body(Z1, N1, Z2, N2, n):

    A1 = Z1 + N1
    A2 = Z2 + N2

    A = A1 + A2 + n
    N = N1 + N2 + n
    Z = Z1 + Z2

    return -(bethe(Z, N) - (bethe(Z1,N1) + bethe(Z2, N2)))

def E2Body_AME(Z1, N1, Z2, N2, n):

    A1 = Z1 + N1
    A2 = Z2 + N2

    A = A1 + A2 + n
    N = N1 + N2 + n
    Z = Z1 + Z2

    return -(binding(Z, N, nominal = 1)*(Z+N) - (binding(Z1,N1, nominal = 1)*(Z1+N1) + binding(Z2, N2, nominal = 1)*(Z2+N2)))/1e3


def E3Body(Z1, N1, Z2, N2, Z3, N3):
    A1 = Z1 + N1
    A2 = Z2 + N2
    A3 = Z3 + N3

    A = A1 + A2 + A3
    N = N1 + N2 + N3
    Z = Z1 + Z2 + Z3

    return bethe(Z, N) - (bethe(Z1,N1) + bethe(Z2, N2) + bethe(Z3, N3))


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from ame import *
from plotinit import initMPL
#initMPL()

ZZ, NN = np.mgrid[1:92+1, 1:144+1]
AA = ZZ + NN


try:
    Ys = np.load("yields.npy")
except:
    Ys = yields(ZZ,NN)
    np.save("yields.npy", Ys)

#E_ame /= 1000 # ame data in MeV
#E_ame[E_ame < 1e-2] = np.nan






# from scipy.optimize import minimize
# E2 = lambda x: -E2Body(x[0], x[1], 92-x[0], (236-92)-x[1]-x[2], x[2])
# bnds = ((0, None), (0, None), (0, 4))
# res = minimize(E2, [236/4, 236/4, 2], method='L-BFGS-B', bounds = bnds, tol=1e-6)
#


# use real data




def EgainA(A1, n = 0, AME = False):
    Z1 = 92/2   # need to define
    N1 = A1-Z1

    A2 = A-A1-n
    Z2 = Z-Z1
    N2 = N-N1-n
    if AME:
        return E2Body_AME(Z1, N1, Z2, N2, n = n)
    else:
        return E2Body(Z1, N1, Z2, N2, n = n)

def EgainZ(Z1, n = 0, AME = False):   # need to define
    N1 = (236-92)/2
    A1 = Z1+N1

    A2 = A-A1-n
    Z2 = Z-Z1
    N2 = N-N1-n
    if AME:
        return E2Body_AME(Z1, N1, Z2, N2, n = n)
    else:
        return E2Body(Z1, N1, Z2, N2, n = n)

def EgainZN(Z1, N1, n = 0, AME = False):   # need to define
    A1 = Z1+N1
    # if A1 > A or N1 > N or Z1 > Z: return np.nan

    A2 = A-A1-n
    Z2 = Z-Z1
    N2 = N-N1-n

    if AME:
        return E2Body_AME(Z1, N1, Z2, N2, n = n)
    else:
        return E2Body(Z1, N1, Z2, N2, n = n)


# Z1 = np.arange(236)[30:63]
# plt.plot(Z1, EgainZ(Z1, n = 0), label = "no neutron")
# plt.plot(Z1, EgainZ(Z1, n = 1), label = "1 neutron")
# plt.plot(Z1, EgainZ(Z1, n = 2), label = "2 neutrons")
#
# Z11 = np.arange(236)[43:50]
# plt.plot(Z11, EgainZ(Z11, n = 0, AME = 1), lw = 3, c = "tab:blue", alpha = 0.5)
# plt.plot(Z11, EgainZ(Z11, n = 1, AME = 1), lw = 3, c = "tab:orange", alpha = 0.5)
# plt.plot(Z11, EgainZ(Z11, n = 2, AME = 1), lw = 3, c = "tab:green", alpha = 0.5)
#
# plt.legend()
# plt.xlabel("proton number $Z_1$ of first fragment")
# plt.ylabel("energy release in the process in MeV")
# plt.axhline(0, ls = "--", color = "gray", zorder = -1)
# plt.axvline(92/2, ls = "--", color = "gray", zorder = -1)
# plt.savefig("bindingBethe.pdf")
f = 7
fig = plt.figure(figsize = (42/f,27/f))
ax = plt.gca()

# E = Ys*0 + np.nan
#
# for z in range(Ys.shape[0]):
#     for n in range(Ys.shape[1]):
#         if ZZ[z,n] < Z and NN[z,n] < N:
#             E[z,n] = EgainZN(ZZ[z,n],NN[z,n])
# plt.matshow(E)

print("energy Cs:", EgainZN(52, 134-55, n = 2, AME = 1))
print("energy Cs + Q0:", EgainZN(52, 134-55, n = 2, AME = 1) + 6.545)

print("energy symm:", EgainZN(55, 134-55, n = 2, AME = 1))
print("energy symm + Q0", EgainZN(46, 118-46, n = 0, AME = 1) + 6.545)

print("energy nucl:", EgainZN(56, 141-56, n = 2, AME = 1))
print("energy nucl + Q0", EgainZN(56, 141-56, n = 2, AME = 1) + 6.545)

import matplotlib as mpl
cmap = mpl.cm.get_cmap("Oranges")


for z in range(Ys.shape[0]):
    for n in range(Ys.shape[1]):

        if Ys[z,n] is not np.nan and Ys[z,n] > 1e-2 and ZZ[z,n] < Z and NN[z,n] < N:

            nn = 2
            print(ZZ[z,n], NN[z,n],EgainZN(ZZ[z,n],NN[z,n], n = nn))


            plt.scatter([NN[z,n]], [ZZ[z,n]], s = Ys[z,n]**2*20000, color = cmap(np.abs(EgainZN(ZZ[z,n],NN[z,n], n = nn, AME = 1)-158)/50))
            #plt.scatter([ZZ[z,n]], [NN[z,n]], s = -EgainZN(ZZ[z,n],NN[z,n])/1000)


plt.scatter(141-56, 56, marker = "o", facecolors=cmap(np.abs(EgainZN(56,141-56, n = nn, AME = 1)-158)/50), edgecolors = "black")
plt.scatter(118-46, 46, marker = "o", facecolors=cmap(np.abs(EgainZN(46,118-46, n = 0, AME = 1)-158)/50), edgecolors = "tab:blue")
plt.ylabel("proton number $Z_1$")
plt.xlabel("neutron number $N_1$")

xl = ax.get_xlim()
yl = ax.get_ylim()

major_ticksx = np.arange(xl[0], xl[1], 5)
minor_ticksx = np.arange(np.floor(xl[0]), np.ceil(xl[1]), 1)

major_ticksy = np.arange(yl[0], yl[1], 5)
minor_ticksy = np.arange(np.floor(yl[0]), np.ceil(yl[1]), 1)

#ax.set_xticks(major_ticksx)
ax.set_xticks(minor_ticksx, minor=True)
#ax.set_yticks(major_ticksy)
ax.set_yticks(minor_ticksy, minor=True)
ax.grid(which='minor', alpha = 0.2, zorder = -1000)

#ax.axis('equal')

#ax.set_xlim(xl)
#ax.set_ylim(yl)

plt.savefig("yields.pdf")
