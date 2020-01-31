# ex 1.1
from uncertainties import ufloat
from ame import excess, binding, mass

mn = ufloat(8071.31713,0.00046)
mn = excess(Z = 0, N = 1, nominal = False)
print(mn)

mp = ufloat(7288.97061,0.00009)
mp = excess(Z = 1, N = 0, nominal = False)

ma = ufloat(2424.91561,0.00006)
ma = excess(Z = 2, N = 2, nominal = False)

dm = 2*mp + 2*mn - ma
print("1.1", dm)
print("1.1", binding(Z=2, N=2, nominal = False)*4)

# ex 1.2

mu238 = ufloat(47307.783, 1.493)
mu238 = excess(Z = 92, N = 238 - 92, nominal = False)

mu234 = ufloat(42444.644, 1.113)
mu234 = excess(Z = 90, N = 234 - 90, nominal = False)

dm = mu238 - (mu234 + ma)
print("1.2", dm)

# TODO uncertainties
#print(mass(90, 236-90))
Ma = ufloat(3.727379378e9, 3.000000023)/1e6 # keV
Ma = mass(Z = 2, N = 2, nominal = False)
Mu236 = ufloat(219.87505e9, 1.86)/1e6 #todo # keV
Mu236 = mass(Z = 92, N = 236-92, nominal = False)
Q = dm


""" energy calculation """

"""
def f(p):
    Ma = 3.727e3
    Mu236    = 219.875e3
    Q = 2438.2
    return Q + Ma + Mu236 - (p**2+Ma**2)**0.5 - (p**2+Mu236**2)**0.5


from scipy.optimize import fsolve
pa = fsolve(f, x0 = 1)[0]

"""
mm2 = Ma**2 + Mu236**2 - (Ma + Mu236 + Q)**2
pa2 = 1/4 * (mm2**2 - 4*Ma**2*Mu236**2) / (Q + Ma + Mu236)**2

print(pa2**0.5)

pa = pa2**0.5
Ekin_e = (pa**2 + Ma**2)**0.5 - Ma
print("exact ", Ekin_e)


""" comparison with simple formula """
Ekin_a = Q * (Mu236)/(Mu236 + Ma)
print("approximation ", Ekin_a)

print("The error made is ", (Ekin_a-Ekin_e)/ ((Ekin_a+Ekin_e)/2))


""" 1.3 """

Elast = mass(92, 235-92, nominal = False) + mass(0,1, nominal = False) - mass(92, 236-92, nominal = False)
print("Binding energy of last neutron in U236: ", Elast)
