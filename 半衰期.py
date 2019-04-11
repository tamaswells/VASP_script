from __future__ import division
from math import exp,log
from scipy import constants as C
Tem=273.15+500
h_p =  6.62606957E-34 # J*s Plank Constant
k_b =  1.38064852E-23 # m²*kg*s⁻²*K^-1 Boltzman Constant 
R_gas =  8.3144598      # J*mol^-1*K^-1 Gas Constant 
l_s =  299792458      # light speed m * s ^-1
kB = C.physical_constants["Boltzmann constant in eV/K"][0]
h = C.physical_constants["Planck constant in eV s"][0]
h_p=h
k_b=kB
Tem = float(Tem)  # Temperature 
beta = 1/(k_b * Tem)
G=70*27.2116/627.51  #kcal/mol->eV
k=1/beta/h_p*exp(-1.0*G*beta)
t_0_5=-log(0.5)/k
print("半衰期为%5.3fh" %(t_0_5/3600))