# define networks

""" NETWORK 1 """

import numpy as np

dim = 5
L = 1

# coupling
J1 = np.zeros((dim, dim))
J1[1,0] = -1 # influence from protein 0 on gene 1
J1[4,0] = +1
J1[3,0] = +1

J1[2,1] = +1

J1[3,2] = -1
J1[4,2] = -1

J1[3,3] = +1
J1[4,3] = +1

J1[4,4] = -1
J1[3,4] = -1

# diffusion
D1 = np.zeros(dim)
D1[4] = 0.4




""" NETWORK 2 """

dim = 5

# coupling
J2 = np.zeros((dim, dim))
J2[1,0] = +1 # influence from protein 0 on gene 1

J2[3,1] = -1
J2[2,1] = +1
J2[4,1] = +1

J2[2,2] = -1
J2[3,2] = +1
J2[4,2] = -1

J2[4,3] = -1

J2[4,4] = +1
J2[3,4] = +1

# diffusion
D2 = np.zeros(dim)
D2[3] = 0.4




""" NETWORK 3 """

dim = 5

# coupling
J3 = np.zeros((dim, dim))
J3[3,0] = +1 # influence from protein 0 on gene 1

J3[4,1] = +1

J3[1,2] = -1
J3[3,2] = -1
J3[4,2] = -1

J3[2,3] = -1
J3[4,3] = +1

J3[4,4] = -1
J3[3,4] = -1
J3[2,4] = +1

# diffusion
D3 = np.zeros(dim)
D3[4] = 0.4


""" NETWORK 4 """

dim = 5

# coupling
J4 = np.zeros((dim, dim))
J4[1,0] = -1 # influence from protein 0 on gene 1

J4[0,1] = +1
J4[2,1] = -1
J4[4,1] = -1

J4[1,2] = +1
J4[3,2] = -1

J4[2,3] = +1
J4[3,3] = +1
J4[4,3] = -1

J4[3,4] = +1

# diffusion
D4 = np.zeros(dim)
D4[2] = 0.4


""" NETWORK R: Repressilator """

# coupling
dim = 3
JR = np.zeros((dim, dim))

JR[1,0] = -1
JR[2,1] = -1
JR[0,2] = -1

# diffusion
DR = np.zeros(dim)


""" NETWORK T: Turing loop """

# coupling
dim = 2
JT = np.zeros((dim, dim))

JT[1,0] = +1
#JT[0,0] = +1

JT[0,1] = -1
#JT[1,1] = -1

# diffusion
DT = np.zeros(dim)
DT[1] = 1 # set the inhibitor to be diffusive


""" hierarchical network H"""
JH = np.zeros((12,12))
JH[2,0] = -1
JH[5,0] = +1
JH[9,0] = +1

JH[3,2] = -1
JH[5,2] = -1
JH[4,2] = -1
JH[9,2] = -1

JH[3,3] = -1
JH[0,3] = -1

JH[0,4] = -1
JH[2,4] = +1

JH[5,5] = -1
JH[7,5] = -1

JH[5,6] = +1

JH[7,7] = +1
JH[8,7] = +1
JH[6,7] = +1

JH[11,8] = -1

JH[9,9] = -1
JH[11,9] = -1

JH[9,10] = +1

JH[10,11] = +1
JH[11,11] = +1

DH = np.zeros(12)
DH[5] = 0.4
DH[6] = 0.4
DH[10] = 0.4
DH[9] = 0.4
