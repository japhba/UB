import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np


def resonance():
    plt.rc('text.latex', preamble=r'\usepackage{nicefrac}')

    x = np.linspace(0,1, 100)
    y = np.sin(x)

    fig = plt.figure(figsize = (10,3))

    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)


    lmbd = 1 # wavelength
    Lmax = 2*lmbd

    ax2.set_ylim(0, Lmax/lmbd)
    ax2.set_xlim(0,200)
    ax2.set_ylabel("size of cell in $\lambda$")
    ax2.set_xlabel("elapsed time $t$")

    ax2.set_yticks([0,0.5,1,1.5,2])
    ax2.set_yticklabels([r"$0$", r"$\frac{1}{2}$", r"$1$", r"$\frac{3}{2}$", r"$2$"])

    ax1.set_xticks([0,0.5,1,1.5,2])
    ax1.set_xticklabels([r"$0$", r"$\frac{1}{2}$", r"$1$", r"$\frac{3}{2}$", r"$2$"])

    ax1.set_ylabel("amplitude of oscillation")
    ax1.set_xlabel("size of cell in $\lambda$")

    plt.tight_layout()

    lWave,   = ax1.plot([0,Lmax],[-1,1])
    #rWave,   = ax1.plot([0,Lmax],[-1,1])
    sWave,   = ax1.plot([0,Lmax],[-1,1])
    #averg,   = ax1.plot([0,2*np.pi],[-1,1])


    T = 1000

    LLine, = ax2.plot([0,T], [0,0])

    k = 2*np.pi / lmbd

    # frequency
    w = 1

    L = 1 # length of resonator

    def sin(x, t):
        return np.sin(x-t)

    def animate(t):

        xx = x * np.max(Lmax) / T * (t+1) # cell growth

        tscale = 10
        t = t/tscale
        c = w/k
        L = np.max(xx)

        y = xx*0

        nRefl = 10

        for j in range(nRefl):
            y += 0.5*np.sin(k*xx + 2*j*k*L)
            y += 0.5*np.sin(-k*(L-xx) + (2*j+1)*k*L)

        lWave.set_data(xx, y/nRefl)
        #rWave.set_data(xx, np.sin(-k*xx-w*(t-(L-xx)/c)))
        sWave.set_data(xx, (y/nRefl)**2)
        #averg.set_data(xx, 0.5*(np.cos(2*k*xx+2*w*(xx-L)/c)+1))

        LLine.set_data([0,t], [0,np.max(xx)])


    ani = matplotlib.animation.FuncAnimation(fig, animate, frames=T, interval = 1000/30, blit = False)
    return ani
