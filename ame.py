import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from uncertainties import ufloat

def bethe(Z,N):

    """ returns total binding energy in MeV """

    A = N+Z

    #if A == 0: return np.nan

    # values based on Benzaid et al., 2019
    aV = 14.9297
    aS = 15.0580
    aC = 0.6615
    aA = 21.6091
    aP = 10.1744

    Av = aV*A
    As = -aS*A**(2/3)
    Ac = -aC*(Z-0)*Z / A**(1/3)
    Aa = -aA*(A-2*Z)**2 / A



    Ap  = np.where(np.logical_and(N % 2 == 0,  Z % 2 == 0), +aP/A**(1/2), 0)
    Ap += np.where(np.logical_and(N % 2 == 1,  Z % 2 == 1), -aP/A**(1/2), 0)

    """
    if N % 2 == 0 and Z % 2 == 0:
        Ap = +aP

    if N % 2 == 1 and Z % 2 == 0:
        Ap = 0

    if N % 2 == 0 and Z % 2 == 1:
        Ap = 0

    if N % 2 == 1 and Z % 2 == 1:
        Ap = -aP
    """

    return Av + As + Ac + Aa + Ap


def excess(Z,N, nominal = False):
    df = pd.read_csv('/Users/jan/Google Drive/Uni/Uebungen/7. Semester/Nuclear Physics/ame2016.csv')
    row = df.loc[(df['Z'] == Z) & (df['N'] == N)]
    try:
        if nominal:
            return ufloat(row["massexcess"].to_numpy()[0], row["uncmassex"].to_numpy()[0]).n
        else:
            return ufloat(row["massexcess"].to_numpy()[0], row["uncmassex"].to_numpy()[0])
    except:
        return np.NaN

def binding(Z,N, nominal = False):
    df = pd.read_csv('/Users/jan/Google Drive/Uni/Uebungen/7. Semester/Nuclear Physics/ame2016.csv')
    row = df.loc[(df['Z'] == Z) & (df['N'] == N)]
    try:
        if nominal:
            return ufloat(row["binding"].to_numpy()[0], row["uncbind"].to_numpy()[0]).n
        else:
            return ufloat(row["binding"].to_numpy()[0], row["uncbind"].to_numpy()[0])
    except:
        return np.NaN


def mass(Z,N, nominal = False):
    return excess(Z,N, nominal) + (Z+N)*931494.102 # atomic mass unit u in keV

excess = np.vectorize(excess, signature = '(),(),()->()', excluded = 'nominal')
binding = np.vectorize(binding, signature = '(),(),()->()', excluded = 'nominal')
mass = np.vectorize(mass, signature = '(),(),()->()', excluded = 'nominal')
