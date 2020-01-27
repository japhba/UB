### CELL MODEL ###

"""

x00: lowercase letters for initial conditions
x0:  lowercase letters for environment variables
x:   lowercase letters for concentration variables

"""

import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
import matplotlib.animation as animation


def plotGraphs():

    t = np.linspace(0, 100, 1000)
    t = [0, 10]

    # We assume that the cell is in a thermodynamic bath, whose concentration of different substances is unaffected by the cell's intake.
    #V = n
    #n = 1/V
    d = 3

    s0 = 1

    V00 = 1
    s00 = s0
    n00 = 0

    x00 = [V00, s00, n00]

    gD = 1
    gn = 1
    gp = 1
    gphi = .1

    phi0 = 1

    from scipy.integrate import solve_ivp

    def get_dx(t,x):
        V = x[0]
        s = x[1]
        n = x[2]

        dV = gn*n
        ds = gn*n + gp*(s0-s) - 1*s
        dn = gphi*phi0 - 1*n

        return [dV, ds, dn]


    sol = solve_ivp(get_dx, (0,100), x00)


    fig = plt.figure(figsize = (10,3))
    gs = fig.add_gridspec(1, 3)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[0,2])

    axes = [ax1, ax2, ax3]
    ttl = ["$V$", "$s$", "$n$"]

    for ax, xx, title in zip(axes, sol.y, ttl):
        ax.plot(sol.t, xx)
        ax.set_title(title)

    plt.show()
