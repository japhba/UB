""" make a calculation for the ratio """

import numpy as np

th = np.arcsin(0.314**0.5)
sigma = 1e-48
dm12 = 7.5e-5 # KAMLAND
dm12 = 4.9e-5
dm123 = 2.5e-3
Ev = 5e6

hc = 3.126e-26
e  = 1.602e-19

meter = e/hc

N = 3.5e51 / 2 # particles in earth
NA = 6.022e23

L_ = 12742e3
V_ = 1.08e21

L = 12742e3*meter
V = 1.08e21*meter**3

n = N / V

G = 1.43585e-62 / (hc)**3 * (e)**2


phi = 2.35e6*10**4
MK = 50000e3
u = 1.66e-27
mH20 = 18.015*u

nsl =  N/V_*sigma*L_

nK = MK/mH20*9
print("Particles in SuperK: ", nK)
print("Particle density in earth: ", N/V_)
print("nSigmaL ", nsl)

print("daynight earth abs", nK*nsl*sigma*phi)
print("daynight earth rel", 2*(1.0-np.exp(-nsl*sigma*phi))/(1+np.exp(-nsl*sigma*phi)))

lv = 4*np.pi*Ev/dm123
l0 = 2**0.5*np.pi/(G*n)

thm = np.arctan(np.sin(2*th) / (np.cos(2*th) + lv/l0))

lm = lv/(1+2*(lv/l0)*np.cos(2*th)+(lv/l0)**2)

"""
alternative solution with simpler formalism
"""

a = 2*2**0.5*G*n*Ev/dm123

C = ((np.cos(2*th)-a)**2+np.sin(2*th)**2)**0.5
Dm2 = dm123*C
sth2 = np.sin(2*th)/C

eps1 = 2**0.5*G*n/(dm12/(2*Ev))


# approximation
# rho in mol / cm**3
rho = N/NA / (V_*1e6)
print("rho ", '{:.2e}'.format(rho))

# big DeltaM

eps2 = rho/100 * 8e-5/dm12 * Ev/5e6

print("eps1 ", '{:.2e}'.format(eps1))
print("eps2 ", '{:.2e}'.format(eps2))

eps = eps1

xi = ((1+eps)**2 - 4*eps*np.cos(th)**2)**0.5


Nd = np.sin(th)**2

Nn = np.sin(th)**2 + 1.7e-2 * (7.37e-5/dm12) * (rho/2) * (Ev/5e6) * np.sin(dm12*L/(4*Ev))**2
print("Nn1", Nn)
#xi = 1
Nn = np.sin(th)**2 * (1 + 4*eps*np.cos(th)**2 * np.sin(xi*dm12*L/(4*Ev))**2 / (xi**2) )
print("Nn2", Nn)

A = 2 * (Nn - Nd) / (Nd + Nn)

print("G ", '{:.2e}'.format(G))
print("lv/l0 ", '{:.2e}'.format(lv/l0))
print("th in °",  '{:.2e}'.format(th/2/np.pi*360))
print("L ",  '{:.2e}'.format(L))
print("ne ",  '{:.2e}'.format(n))
print("V ",  '{:.2e}'.format(V))
print("thm ", '{:.2e}'.format(thm))

print("eps ", '{:.2e}'.format(eps))
print("xi ", '{:.2e}'.format(xi))
print("---")
print("A ", '{:.2e}'.format(A))
