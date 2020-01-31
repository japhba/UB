
from molmass import ELEMENTS
import pandas as pd
import numpy as np

import re

df = pd.read_csv('/Users/jan/Google Drive/Uni/Uebungen/7. Semester/Nuclear Physics/yield.csv', sep=';')
s = df['Product'].to_numpy()

ZS = [re.split('(\d+)',s[i])[0] for i in range(len(s))]
AS = [re.split('(\d+)',s[i])[1] for i in range(len(s))]

zs = [ELEMENTS[ZS[i]].number for i in range(len(s))]
as_ = np.array(AS).astype(int)

ns = as_ - zs

df['Z'] = pd.Series(zs, index=df.index)
df['N'] = pd.Series(ns, index=df.index)
df.to_csv("yields2.csv", sep=';')
