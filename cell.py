import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

gm = 6
gp = 1

#p = np.zeros(dim) # protein expression
#m = np.zeros(dim) # mRNA readout

#D = np.zeros(dim) # diffusion
#T = np.zeros(dim) # threshold

T0 = np.array([-0.01, -0.03, 0.02, 0.01, -0.02]) #default treshhold


class Cell():

    def __init__(self, J, D, m, p, T = T0, name = ""):
        J = np.array(J)
        D = np.array(D)

        self.J = J
        self.D0 = D
        self.T0 = np.array(T)

        self.maxL = 2**5

        self.L = 1 #initial number of cells

        self.m = m
        self.p = p

        self.dim = J.shape[0]

        self.ps = []
        self.ms = []

        self.gen = 0

        self.name = name

        # use the same parameters over all cells
        self.D = np.repeat(self.D0[:, np.newaxis], self.L, axis=1)
        self.T = np.repeat(self.T0[:, np.newaxis], self.L, axis=1)

    def evolve(self, t):
        """
        evolves the cells forward a time t and returns the evolution history over this
        period. The cell's state is then set to the last stage in the evolution.
        """

        p = self.p
        m = self.m

        L = self.L
        dim = self.dim

        D = self.D
        J = self.J
        T = self.T

        def f(x, beta = 40):
            return 1/(1 + np.exp(-beta*x))

        def get_dx(x , times): # x = i x l (flattened)
            x = x.reshape((2*dim, L))

            #if L == 4: print(x)

            p = x[:dim,:]
            m = x[dim:,:]

            P = p.mean(axis = 1)

            # use the same parameters over all cells
            P = np.repeat(P[:, np.newaxis], L, axis=1)


            ### GOVERNING EQUATIONS

            dpdt = gp*(m-p) + D*(P-p)

            F = f( np.einsum('ij,jl->il', J, p) - T )
            dmdt = gm*(F-m)

            ###

            dxdt = np.concatenate((dpdt, dmdt), axis = 0)

            dxdt = dxdt.reshape((2*dim*L))

            return dxdt # is not explicitely time dependent


        times = t
        x_ini = np.concatenate((self.p, self.m), axis = 0).reshape(2*dim*L)

        x = odeint(get_dx, x_ini, times)
        x = x.reshape((len(times), 2*dim, L))

        ps, ms = x[:,:dim,:], x[:,dim:,:]

        # set current state
        self.p = x[-1,:dim,:]
        self.m = x[-1,dim:,:]

        # extend cell time series
        self.ps.append(x[:,:dim,:])
        self.ms.append(x[:,dim:,:])

        return self


    def divide(self, sigma = 1e-5):
        oldL = self.L
        newL = self.L * 2

        dim = self.dim

        # introduce noise
        etaP = np.random.uniform(-sigma, sigma, (dim,oldL))
        etaM = np.random.uniform(-sigma, sigma, (dim,oldL))

        #if newL == 4: print(etaP)


        p1 = self.p + etaP[:, :]
        p2 = self.p - etaP[:, :]

        m1 = self.m + etaM[:, :]
        m2 = self.m - etaM[:, :]

        self.p = np.concatenate((p1, p2), axis = 1)
        self.m = np.concatenate((m1, m2), axis = 1)

        # use the same parameters over all cells
        self.D = np.repeat(self.D0[:, np.newaxis], newL, axis=1)
        self.T = np.repeat(self.T0[:, np.newaxis], newL, axis=1)

        self.L = newL
        self.gen += 1



    def color(self, a=[0.5], j = 1):
        """constructs appropriate color map"""
        gen = self.gen + 1
        if j < gen:
            a2 = np.concatenate((a, a))
            for i in range(len(a)):
                a2[i] = a[i] - 2**-(j+1)
                a2[i+2**(j-1)] = a[i] + 2**-(j+1)
            return self.color(a2, j+1)

        if j == gen:
            return np.array(a)*0.9 + 0.1
