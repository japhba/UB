import numpy as np

""" all energies in MeV """

# single particle energies
E1 = -4.15
E2 = -3.28
E3 = +0.93

# interaction
V11 = -2.0094
V12 = -1.3225
V13 = -3.8935
V22 = -2.3068
V23 = -0.8385
V33 = -0.8119

V21 = +V12
V32 = +V23
V31 = +V13

E = np.array([[E1,0,0],[0,E2,0],[0,0,E3]])
V = np.array([[V11,V12,V13],[V21,V22,V23],[V31,V32,V33]])

H = 2*E + V

# compare with experimental data
ev, U = np.linalg.eig(H)

from ame import excess, mass
experiment = excess(Z = 8, N = 18 - 8, nominal = 0) - (excess(Z = 8, N = 16 - 8, nominal = 0) + 2*excess(Z = 0, N = 1, nominal = False))

i = np.argmin(ev)

print("Ground state: ", [round(ev[i], 2) for i in range(3)]) # MeV
print("Ground state: ", round(ev[i], 2)) # MeV
print("Measurement: ", experiment/1e3) # from keV to MeV

print("transition elements:")
print("<1|GS>:", (U[0,i]))
print("<2|GS>:", (U[1,i]))
print("<3|GS>:", (U[2,i]))

print("probabilities:")
print("<GS|1>:", np.abs(U[0,i])**2)
print("<GS|2>:", np.abs(U[1,i])**2)
print("<GS|3>:", np.abs(U[2,i])**2)
