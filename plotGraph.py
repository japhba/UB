import matplotlib.pyplot as plt
import networkx as nx                                 #ネットワーク用
import matplotlib.animation as animation              #アニメーション用
import matplotlib.patches as patches                  #円とかの描画用
from matplotlib.patches import FancyArrowPatch

from scipy.integrate import solve_bvp, odeint

from ipywidgets import FloatSlider

import numpy as np

import matplotlib
# get 5 different colormaps for each protein
cm1 = matplotlib.cm.get_cmap('Purples')
cm2 = matplotlib.cm.get_cmap('Greens')
cm3 = matplotlib.cm.get_cmap('Oranges')
cm4 = matplotlib.cm.get_cmap('Reds')
cm5 = matplotlib.cm.get_cmap('Blues')
cm6 = matplotlib.cm.get_cmap('PuRd')
cm7 = matplotlib.cm.get_cmap('GnBu')
cm8 = matplotlib.cm.get_cmap('OrRd')

cms = [cm1, cm2, cm3, cm4, cm5, cm6, cm7, cm8]


def plotDevelopment(cell, divs = 3, tmax = 120, proteins = [0], sigma = 1e-5):

    fig = plt.figure(figsize = (10,3))
    gs = fig.add_gridspec(1, 10)
    ax1 = fig.add_subplot(gs[:,:3])
    ax2 = fig.add_subplot(gs[:,4:9])


    dt=0.1

    if type(divs) is int:
        divs += 1
        divs = [i*(tmax // divs) for i in range(divs)]

    if divs[0] is not 0: divs = np.concatenate(([0],divs))

    divs = np.concatenate((divs, [tmax]))

    cell.divs = divs


    for d in range(len(divs)-1): # plot every generation

        # evolve net
        t=np.arange(divs[d],divs[d+1],dt)
        cell.evolve(t-divs[d])

        # plot epoch
        for l in range(cell.gen+1):
            for i in proteins:
                ax2.axvline(divs[d+1], zorder = -1, c = "xkcd:light gray")
                label = '$p_'+str(i)+'$' if l is 0 and d is 0 else None
                ax2.plot(t, cell.ps[cell.gen][:,i,l], c = cms[i](cell.color()[l], alpha = 0.8), lw = 2/2**d, label = label)
                ax2.set_xlabel("time")
                ax2.set_ylabel("protein abundance")

                xT=ax2.get_xticks()
                xL=ax2.get_xticklabels()

        # divide cell
        cell.divide(sigma = sigma)

    legend = plt.figlegend(loc='center right')
    title  = ax2.set_title(str(cell.name))

    from plotGraph import animateGraph
    ani = animateGraph(fig, ax1, cell, times = np.arange(0, tmax, dt), ax2 = ax2, xL = None) # weird behavior with x ticks

    plt.tight_layout()

    return ani, cell
    #plt.close()

def plotProj(cell):
    plt.close()
    fig = plt.figure(figsize = (4,4))
    from sklearn.manifold import TSNE
    from sklearn.decomposition import PCA
    X = np.concatenate(tuple(np.concatenate(tuple(cell.ps[g][:,:,l] for l in range(cell.ps[g].shape[2]))) for g in range(len(cell.ps))))
    #X_embedded = TSNE(n_components=2).fit_transform(X)
    X_embedded = PCA(n_components=2).fit_transform(X)

    gens = len(cell.ps)
    ms = ["o", "o", "o"]

    colors = np.concatenate(tuple(np.concatenate(tuple(cms[l % 7](np.linspace(g/gens, (g+1)/gens, len(cell.ps[g][:,:,l]))) for l in range(cell.ps[g].shape[2]))) for g in range(gens)))
    #marker = np.concatenate(tuple(np.concatenate(tuple([ms[l] for _ in range(len(cell.ps[g][:,:,l]))] for l in range(cell.ps[g].shape[2]))) for g in range(gens)))
    plt.scatter(X_embedded[:,0], X_embedded[:,1], c = colors, alpha = .1, s = 2)
    plt.title("Phase space orbit")
    plt.xlabel("PCA component 1")
    plt.ylabel("PCA component 2")

    plt.tight_layout()

    #plt.gca().set_aspect('equal', 'box')



def animateGraph(fig, ax, cell, times, ax2 = None, xL = None):

    """ construct plot """

    n = cell.dim
    J = cell.J
    D = cell.D0

    x = np.concatenate(tuple(cell.ps[g][:,:,0] for g in range(len(cell.ps)))) #time series

    Graph = nx.complete_graph(n, create_using=nx.DiGraph())
    pos   = nx.spring_layout(Graph, k=1)

    node_colors = ["tab:purple", "tab:green", "tab:orange", "tab:red", "tab:blue"]
    node_edge_colors = ["black" if D == 0 else "xkcd:light grey" for D in D]

    edge_colors = [['tab:green' if J[i,j] > 0 else 'tab:red' for j in range(n)] for i in range(n)]

    ax.axis('off')

    """ animation """
    skip=10

    t_max = times[-1]
    dt = times[1]-times[0]

    def animate(i, Graph, x, w, pos, skip):

        # option 2, remove all lines and collections
        for artist in ax.lines + ax.collections + ax.texts:
            artist.remove()

        # self feedback
        for node in Graph:
            ax.add_patch(patches.Circle(xy=(pos[node][0]-0.1,pos[node][1]),
                                        radius=0.15, fill=False,
                                        ec=edge_colors[node][node],
                                        linewidth=abs(w[node,node]*x[i*skip,node]*5),
                                        alpha=1, zorder = -2))

            ax.text(pos[node][0]-0.02, pos[node][1]-0.02, node, color = "xkcd:light gray")

        nx.draw_networkx_nodes(Graph, pos, node_color=node_colors[:n],
                               alpha=1,
                               node_size=abs(x[i*skip,:]*1000), # DO VALUES GET NEGATIVE??
                               linewidths=2.0,
                               edgecolors=node_edge_colors,
                               ax = ax)


        for (u,v) in Graph.edges():
            if w[v,u] != 0:
                ax.add_patch(FancyArrowPatch(pos[u],pos[v],arrowstyle='-|>', # tail, head
                                             connectionstyle='arc3,rad=0.1',
                                             mutation_scale=10.0,
                                             lw=abs(w[v,u]*x[i*skip,u]*5),
                                             alpha=1,
                                             color=edge_colors[v][u],
                                             shrinkA=20,shrinkB=10))


        ax.text(0.9, 0.9, "$t =$" + str(int(i*skip*dt)), transform=ax.transAxes)

        if ax2 is not None:
            ax2.tick_params(
                            axis='both',
                            which='both',
                            bottom=True,
                            left=True,
                            labelbottom=True,
                            labelleft=True)

            if xL is not None: ax2.set_xticklabels(xL)

    ani = animation.FuncAnimation(fig, animate, fargs=(Graph,x,J,pos,skip), interval = 400, frames=len(times) // skip, blit = False)
    return ani

def proteinsInCells(cell):
    fig = plt.figure(figsize = (10,10))
    I,J = 4, 4
    gs = fig.add_gridspec(I, J)
    L = 2**(cell.gen-1)

    dt = 0.1

    axes = np.array([[fig.add_subplot(gs[i,j]) if (i)*J + (j+1) <= L else None for j in range(J)] for i in range(I)])

    norms = []
    for l in range(L):
        norms.append((cell.ps[-1][:,:,l]**2).sum(axis=0).sum(axis = 0))


    for ax, l in zip(axes.flatten()[:L], range(L)):
        for i in range(cell.dim):
            t = np.arange(cell.divs[-2], cell.divs[-1], dt)
            ax.plot(t, cell.ps[-1][:,i,l], color = cms[i](0.5), lw = 1, alpha = 0.5)

            if norms[l] < np.mean(norms):
                ax.set_title("Cell " + str(l+1) + " (stem)")

            if norms[l] > np.mean(norms):
                ax.set_title("Cell " + str(l+1) + " (diff)")

    fig.suptitle('Protein abundances in cell ensemble', fontsize=14)

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

cx = None
cy1 = None
cy2 = None

S = None



def phasePlot2D(J = None, D = None):
    fig, ax3 = plt.subplots(figsize=(5,5))
    X, Y = np.meshgrid(np.linspace(0, 1, 20), np.linspace(0, 1, 20))
    x, y = np.meshgrid(np.linspace(0, 1, 500), np.linspace(0,1,500))

    ## apply function
    if J is None:
        J = np.zeros((2,2))
        J[0,0] = 1
        J[1,0] = 1
        J[1,1] = 0
        J[0,1] = -1
    else: J = np.array(J)

    if D is None:
        D = np.array([0, .14])
    else: D = np.array(D)

    def grad(X, Y, I = 1):
        T = np.array([-0.1, 0.2])
        if len(X.shape) > 0:
            T = np.repeat(T[:, np.newaxis], X.shape[0], axis=1)
        if len(X.shape) > 1:
            T = np.repeat(T[:, np.newaxis], X.shape[0], axis=1)

        XY = np.array([X,Y])

        f = lambda x: (1+np.exp(-40*x))**(-1)
        #f = lambda x: np.heaviside(x, 0) # overflow?
        f_arg = np.einsum('ij,j...->i...', J, XY) - T
        if f_arg.any() > 3: print(f_arg)

        F = f( f_arg )

        DX = F[0] - X + I*D[0]*(Y-X)
        DY = F[1] - Y + I*D[1]*(X-Y)

        M = np.hypot(X, Y)

        return DX, DY, M

    DX1, DY1, M1 = grad(X,Y, I = 1)
    dx1, dy1, m1 = grad(x,y, I = 1)

    DX2, DY2, M2 = grad(X,Y, I = -1)
    dx2, dy2, m2 = grad(x,y, I = -1)

    global S
    S = ax3.streamplot(X, Y, DX1, DY1, color = M1) #alpha?


    def limitCycleFinder(): #credit Lutz Lehmann, math.stackexchange, q. 2679760

        def odesys(u,t): return grad(u[0], u[1])[:2]
        def bc(ya, yb): return yb-ya

        def norm(a): return max(abs(a))

        points = [np.array([10.0,10.0])]; # have some non-relevant point

        n = 10
        t_init = np.linspace(0, 2*n*np.pi, n*10+1)
        t_sol  = np.linspace(0, 2*n*np.pi, n*500+1)

        for x in np.linspace(0.2, 0.25, 5+1):
          for y in np.linspace(0.3, 0.35, 4+1):
            u0 = np.array([ x, y ]);
            u_init = odeint(odesys, u0, t_init);
            #print(u_init)

            res = solve_bvp(lambda t,u:odesys(u,t), bc, t_init, u_init.T, tol=1e-2, max_nodes=n*500+1)
            if res.success:
                return res.sol(t_sol)

            #print res.message,"\n",n,": ",res.sol(0)
            if res.success and min( norm(pp - res.sol(0)) for pp in points ) > 5e-3:
                res = solve_bvp(lambda t,u:odesys(u,t), bc, res.x, res.y, tol=1e-3, max_nodes=n*10000+1)
                #print res.message, "\n refine to",res.sol(0)
                if res.success and min( norm(pp - res.sol(0)) for pp in points ) > 1e-2:
                    for j in range(n): points.append(res.sol(2*np.pi*j))
                    u_sol = res.sol(t_sol);

        return u_sol

    """ does not work in a reproducible way

    u_sol = limitCycleFinder()
    ax3.plot(u_sol[0], u_sol[1])

    """

    # nullclines
    global cx, cy1, cy2
    cx = ax3.contour(x, y, dx1, levels=[0], linewidths=2, colors='tab:blue')
    cy1 = ax3.contour(x, y, dy1, levels=[0], linewidths=2, colors='tab:orange')
    #cy2 = ax3.contour(x, y, dy2, levels=[0], linewidths=2, colors='tab:orange')

    ax3.set_xlabel("$x$", color = "tab:blue")
    ax3.set_ylabel("$y$", color = "tab:orange")

    ax3.set_title("Phase space portrait")
    #Q = ax3.quiver(X, Y, DX, DY, M, units='x', pivot='mid')

    #qk = ax3.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',coordinates='figure')
    #ax3.scatter(X, Y, color='0.5', s=1)

    def ISlider(I=1):
        global cx, cy1, cy2, S

        I = kSliderWidget.value
        x, y = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0,1,100))
        dx1, dy1, m1 = grad(x,y, +I)
        dx2, dy2, m2 = grad(x,y, -I)


        # clear lines
        for coll in cx.collections:
            coll.remove()

        for coll in cy1.collections:
            coll.remove()


        #for coll in cy2.collections:
        #    coll.remove()

        S.lines.remove()

        keep = lambda x: not isinstance(x, matplotlib.patches.FancyArrowPatch)
        ax3.patches = [patch for patch in ax3.patches if keep(patch)]

        S = ax3.streamplot(X, Y, DX1, DY1, color = M1) #alpha?



        cx = ax3.contour(x, y, dx1, levels=[0], linewidths=2, colors='tab:blue')
        cy1 = ax3.contour(x, y, dy1, levels=[0], linewidths=2, colors='tab:orange')
        #cy2 = ax3.contour(x, y, dy2, levels=[0], linewidths=2, colors='tab:orange')

        DX1, DY1, M1 = grad(X,Y, I = I)



        fig.canvas.draw_idle()


    kSliderWidget = FloatSlider(
        value=1,
        min=-10.0,
        max=10.0,
        step=0.1,
        continuous_update=False,
        description='$\\alpha$')
    kSliderWidget.observe(ISlider)

    display(kSliderWidget)

    plt.show()

def bifurcation(J = None, D = None):

    fig = plt.figure(figsize = (8,8))
    gs = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])

    axes = [ax1, ax2, ax3, ax4]

    ## apply function
    if J is None:
        J = np.zeros((2,2))
        J[0,0] = 1
        J[1,0] = 1
        J[1,1] = 0
        J[0,1] = -1
    else: J = np.array(J)

    def grad(x, d = .14):
        T = np.array([[-0.1, -0.1], [0.2, 0.2]])

        x1 = x[0]
        y1 = x[1]

        x2 = x[2]
        y2 = x[3]

        X = np.array([[x1, x2], [y1, y2]])

        f = lambda x: (1 + np.exp(-40*x))**(-1)
        #f = lambda x: np.heaviside(x, 0) # overflow?
        f_arg = np.einsum('ij,jl...->il...', J, X) - T

        F = f( f_arg )

        dX = F - X
        dX[1,0] += d*(y2-y1) # save some tensor operations
        dX[1,1] -= d*(y2-y1)

        dx1 = f(x1-y1-(-0.1)) - x1
        dy1 = f(x1-(0.2)) - y1 + d*((y1+y2)/2-y1)

        dx2 = f(x2-y2-(-0.1)) - x2
        dy2 = f(x2-(0.2)) - y2 + d*((y1+y2)/2-y2)

        return [dx1, dy1, dx2, dy2]

    d = np.linspace(0.07, 0.14, 50)
    d_disc = [0.0, 0.12, 0.1225, 0.14]
    xinit = np.random.uniform(0, 1, (2))
    noise = np.random.uniform(0, 1e-3,(2))

    xinit1 = xinit + noise
    xinit2 = xinit - noise

    xinit = np.concatenate((xinit1, xinit2))

    dt = .1
    t = np.arange(0,120,dt)
    sols_d = np.array([odeint(lambda x,t: grad(x, d = dd), xinit, t)  for dd in d_disc])

    for ax, d, sol in zip(axes, d_disc, sols_d):
        # plot the xi
        ax.plot(t, sol[:,0])
        ax.plot(t, sol[:,2])
        ax.set_title("$D = " + str(d) + "$")
        ax.set_xlabel("time")
        ax.set_ylabel("protein expression $x_i$")


    fig.suptitle("Coupling dependence of dynamics", fontsize = 14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()


def Ddependence(J = None, D = None, ds = [0.0, 0.12, 0.1225, 0.14], t = np.arange(0,120, 0.1)):

    d_disc = ds

    fig = plt.figure(figsize = (8,8))
    gs  = fig.add_gridspec(2, 2)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])

    axes = [ax1, ax2, ax3, ax4]

    ## asymmetric simplified model
    if J is None:
        J = np.zeros((2,2))
        J[0,0] = 1
        J[1,0] = 1
        J[1,1] = 0
        J[0,1] = -1
    else: J = np.array(J)

    def grad(x, d = .14):
        T = np.array([[-0.1, -0.1], [0.2, 0.2]])

        x1 = x[0]
        y1 = x[1]

        x2 = x[2]
        y2 = x[3]

        X = np.array([[x1, x2], [y1, y2]])

        f = lambda x: (1 + np.exp(-40*x))**(-1)
        #f = lambda x: np.heaviside(x, 0) # overflow?
        f_arg = np.einsum('ij,jl...->il...', J, X) - T

        F = f( f_arg )

        dX = F - X
        dX[1,0] += d*(y2-y1) # save some tensor operations
        dX[1,1] -= d*(y2-y1)

        dx1 = f(x1-y1-(-0.1)) - x1
        dy1 = f(x1-(0.2)) - y1 + d*((y1+y2)/2-y1)
        #dy1 = f(x1-(0.2)) - y1 + d*(y2-y1)

        dx2 = f(x2-y2-(-0.1)) - x2
        dy2 = f(x2-(0.2)) - y2 + d*((y1+y2)/2-y2)
        #dy2 = f(x2-(0.2)) - y2 + d*(y1-y2)

        return [dx1, dy1, dx2, dy2]


    # random post-splitting starting state
    xinit = np.random.uniform(0, 1, (2))
    noise = np.random.uniform(0, 1e-3,(2))

    xinit1 = xinit + noise
    xinit2 = xinit - noise

    xinit = np.concatenate((xinit1, xinit2))


    sols_d = np.array([odeint(lambda x,t: grad(x, d = dd), xinit, t)  for dd in d_disc])

    for ax, d, sol in zip(axes, d_disc, sols_d):
        # plot the xi
        ax.plot(t, sol[:,0], color = "tab:purple")
        ax.plot(t, sol[:,2], color = "tab:green")
        ax.set_title("$D = " + str(d) + "$")
        ax.set_xlabel("time")
        ax.set_ylabel("protein expression $x_i$")


    fig.suptitle("Coupling dependence of dynamics", fontsize = 14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

def bifurcation(J = None, D = None, supportPoints = 50):

    fig = plt.figure(figsize = (8,4))

    ## apply function
    if J is None:
        J = np.zeros((2,2))
        J[0,0] = 1
        J[1,0] = 1
        J[1,1] = 0
        J[0,1] = -1
    else: J = np.array(J)

    def grad(x, d = .14):
        T = np.array([[-0.1, -0.1], [0.2, 0.2]])

        x1 = x[0]
        y1 = x[1]

        x2 = x[2]
        y2 = x[3]

        X = np.array([[x1, x2], [y1, y2]])

        f = lambda x: (1 + np.exp((-40)*x))**(-1)
        #f = lambda x: np.heaviside(x, 0) # overflow?
        f_arg = np.einsum('ij,jl...->il...', J, X) - T

        F = f( f_arg )

        dX = F - X
        dX[1,0] += d*(y2-y1) # save some tensor operations
        dX[1,1] -= d*(y2-y1)

        dx1 = f(x1-y1-(-0.1)) - x1
        dy1 = f(x1-(0.2)) - y1 + d*((y1+y2)/2-y1)

        dx2 = f(x2-y2-(-0.1)) - x2
        dy2 = f(x2-(0.2)) - y2 + d*((y1+y2)/2-y2)

        return [dx1, dy1, dx2, dy2]

    d = np.linspace(0.07, 0.14, supportPoints)
    xinit = np.random.uniform(0, 1, (2))
    noise = np.random.uniform(0, 1e-3, (2))

    xinit1 = xinit + noise
    xinit2 = xinit - noise

    xinit = np.concatenate((xinit1, xinit2))

    dt = .1
    # late evolution
    #t = np.concatenate(([0],np.arange(100,150,dt)))
    t = np.arange(0,250,dt)
    sols_d = np.array([odeint(lambda x,t: grad(x, d = dd), xinit, t)  for dd in d])

    from scipy.signal import find_peaks
    colors = ["tab:purple", "tab:green"]
    for xx, c in zip([0,2], colors):
        for i in range(len(d)):
            y = sols_d[i, 2000:, xx]
            peaks, _ = find_peaks(y, height=0)
            # if i == 40:
            #     plt.plot(t, y)
            #     plt.plot(t[peaks], y[peaks], ls = "none", marker = "x")
            [plt.scatter([d[i]], [y[peak]], s = 1, color = c, alpha = 0.5) for peak in peaks]

    #max_x1 = np.max(sols_d[:, :, 0], axis = 1)
    #max_x2 = np.max(sols_d[:, :, 2], axis = 1)

    plt.title("Bifurcation in dependence of coupling")
    plt.xlabel("coupling strength $D$")
    plt.ylabel("local maxima of $x_i$ at late times")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
