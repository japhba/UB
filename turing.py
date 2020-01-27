from ipywidgets import *
import matplotlib.pyplot as plt
import matplotlib.animation
import numpy as np

pause = True

def turingAnimation(J, D):

    pause = True

    angle = lambda z: np.angle(z)

    M = lambda k: J + k**2*np.diag(D)

    x = np.linspace(0, 1, 100)
    t = 0

    # figure setup
    fig = plt.figure(figsize = (10,4))
    gs = fig.add_gridspec(10, 3)
    ax1 = fig.add_subplot(gs[2:9, 0:2])
    ax2 = fig.add_subplot(gs[2:9, 2], projection='polar')

    xT=ax2.get_xticks()
    xL=['0',r'$\frac{\pi}{4}$',r'$\frac{\pi}{2}$',r'$\frac{3\pi}{4}$',\
        r'$\pi$',r'$\frac{5\pi}{4}$',r'$\frac{3\pi}{2}$',r'$\frac{7\pi}{4}$']
    ax2.set_xticks(xT)
    ax2.set_xticklabels(xL)

    # initialize second axis
    l1, = ax2.plot([0],[10],marker='o') # (rho, theta)
    l2, = ax2.plot([0],[10],marker='o')

    # solve the Eigenvalue problem for a specific k
    k = 1
    def turingPattern(k=1):
        eValues, eVectors = np.linalg.eig(M(k))
        #print("{num.real:+0.04f} {num.imag:+0.04f}j".format(num=eValues[0]))
        #print("{num.real:+0.04f} {num.imag:+0.04f}j".format(num=eValues[1]))
        l1.set_data([angle(eValues[0])], [np.abs(eValues[0])])
        l2.set_data([angle(eValues[1])], [np.abs(eValues[1])])

        sol = lambda x,t: [[np.real(eVector_i*np.exp(1j*k*x + eValue*t)) for eVector_i in eVector] for eVector, eValue in zip(eVectors, eValues)]
        return sol

    def kSlider(k=1):
        global pause
        #pause = True
        k = kSliderWidget.value
        l = lDropdownWidget.value
        x1.set_ydata(turingPattern(k)(x,t)[l][0])
        x2.set_ydata(turingPattern(k)(x,t)[l][1])
        fig.canvas.draw_idle()

    def lambdaSelector(l = 0):
        global pause
        #pause = True
        k = kSliderWidget.value
        l = lDropdownWidget.value
        x1.set_ydata(turingPattern(k)(x,t)[l][0])
        x2.set_ydata(turingPattern(k)(x,t)[l][1])
        fig.canvas.draw_idle()

    def tEvolution(start):
        global pause # gets the pause variable defined in the outer scope
        pause ^= True

    # initialize the first axis
    x1, = ax1.plot(turingPattern(k=1)(x,t)[0][0], ls = "-")
    x2, = ax1.plot(turingPattern(k=1)(x,t)[0][1], ls = ":")

    ax2.set_title("eigenvalues $\\lambda_{1,2}$", pad = 20)
    ax1.set_title("concentration variation $\\delta c_{1,2}$", pad = 20)
    ax1.set_xlabel("space coordinate $x$")
    ax1.set_ylabel("concentration $\\delta c_{1,2}$")


    # define interactive widgets

    kSliderWidget = FloatSlider(
    value=1,
    min=0,
    max=10.0,
    step=0.1,
    description='$\\kappa$')
    kSliderWidget.observe(kSlider)

    speedSliderWidget = FloatLogSlider(
        value=1e-2,
        base=10,
        min=-10, # max exponent of base
        max=10, # min exponent of base
        step=0.2, # exponent step
    description='speed')

    lDropdownWidget = Dropdown(
        options=[('Lambda 1', 0), ('Lambda 2', 1)],
        value=0,
        description='Eigenmode:',
    )
    lDropdownWidget.observe(lambdaSelector)

    tEvolutionWidget = ToggleButton(
        value=False,
        description='Start Time Evolution',
        disabled=False,
        button_style='', # 'success', 'info', 'warning', 'danger' or ''
        tooltip='Description',
        icon=''
    )
    tEvolutionWidget.observe(tEvolution)

    display(kSliderWidget)
    display(lDropdownWidget)
    display(tEvolutionWidget)
    display(speedSliderWidget)

    # start time animation
    T = 1000
    def animateTuring(t):
        t = t*speedSliderWidget.value
        global pause
        if not pause:
            k = kSliderWidget.value
            l = lDropdownWidget.value
            x1.set_ydata(turingPattern(k)(x,t)[l][0])
            x2.set_ydata(turingPattern(k)(x,t)[l][1])

    ani = matplotlib.animation.FuncAnimation(fig, animateTuring, frames=T, interval = 1000/30, blit = False)

    return ani
