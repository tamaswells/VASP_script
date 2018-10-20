#!/usr/bin/env python27 
# -*- coding: utf-8 -*-
'''
Writen By Qiang Li from ICIQ in 10-03-2018
Modified By Nxu, Jincheng Liu 2018-10-7
This script is used to calculate
1) the vibration contributions to the entropy and enthalpy under specific temperature
2) zero point energies

Only vibrations are considered. No Translation, No Rotation, No Electron contribitions.

Ref: Atkin's Physical chemistry, 10th Edition, Statistical thermodynamics: the vibrational Contribution, Page: 642
	Florian Schlossar Disertation_Appendix_thermo appendix B
To use it: 

python entropy.py 298.15  
 
298.15 is the Temperature, of course you can use other T values in K)

'''
import sys
import math 
#from scipy import constants as con

script, Tem = sys.argv
h_p =  6.62606957E-34 # J*s Plank Constant
k_b =  1.38064852E-23 # m²*kg*s⁻²*K^-1 Boltzman Constant 
R_gas =  8.3144598      # J*mol^-1*K^-1 Gas Constant 
l_s =  299792458      # light speed m * s ^-1
Tem = float(Tem)  # Temperature 
beta = 1/(k_b * Tem)

# get entropy & enthalpy ref: Florian Schlossar Disertation_Appendix_thermo.pdf
def get_pf(miu,index): # get partition function
    x_i = h_p * float(miu) * l_s * beta #βhcv
    if index==0:
        pf_l  = x_i / (math.exp(x_i) - 1)   # Left part in the entropy equation 
        pf_r  = math.log(1 - math.exp(-x_i)) 
        pf    = pf_l - pf_r 
        entropy =  R_gas * pf 
        return entropy
    elif index==1:
        theta_v = x_i * Tem
        enthalpy = R_gas * theta_v * (1.0/(math.exp(theta_v/Tem)-1))
        return enthalpy
    else:
        theta_v = x_i * Tem
        zpes = R_gas * theta_v * 0.5
        return zpes

def collect_miu():
    look_up1 = 'cm-1'
    look_up2 = 'f/i'
    miu_list = []
    zpe_list = []
    with open('OUTCAR') as data_in:
        for i in data_in.readlines():
            if look_up1 in i.rstrip() and look_up2 not in i.rstrip():
                miu_list.append(float(i.split()[7]))
                zpe_list.append(float(i.split()[9]))
    return miu_list, zpe_list

miu_list, zpe_list = collect_miu()

E_zpe = sum(zpe_list)/2000.0  # calculate ZPE 1) Sum(hv)/2  2) convert meV to eV or BY partition function

reduce_l =  [i if i >= 50.0 else 50.0 for i in miu_list] # Convert frequence with wavenum < 50 cm-1 to 50 cm-1,like ASE does !!!!

# Calculate the Entropy S,  i * 100: convert cm-1 to m-1
sum_entropy =  sum(get_pf(i*100,0) for i in reduce_l) 

# convert J * K^-1 * mol^-1 to eV * K^-1 from http://web.utk.edu/~rcompton/constants 
sum_entropy = sum_entropy/1000.0/96.485  

# Calculate the enthalpy H,  i * 100: convert cm-1 to m-1
sum_enthalpy =  sum(get_pf(i*100,1) for i in miu_list) 

# convert J * mol^-1 to eV  from http://web.utk.edu/~rcompton/constants 
sum_enthalpy = sum_enthalpy/1000.0/96.485  # Ev correction to U
KbT = k_b * Tem /(1.60217733e-19) 
#sum_enthalpy+=KbT   # H=U+KbT 

# Calculate ZPE from partition function,  i * 100: convert cm-1 to m-1
sum_ZPE =  sum(get_pf(i*100,2) for i in miu_list) 

# convert J * mol^-1 to eV  from http://web.utk.edu/~rcompton/constants 
sum_ZPE = sum_ZPE/1000.0/96.485 

TS = Tem * sum_entropy  # entropy contribution: T * S 

print "\n			H = E_DFT + E_ZPE + E_H + (nRT  n*%6.7f if in gas phase)\n			G = H - TS = E_DFT + E_ZPE + E_H - T*S\n" %(KbT)
print ' S: %6.7f eV * K^-1   \tE_H: %6.4f eV    \tT*S: %6.4f eV    \tE_ZPE: %6.4f eV' %(sum_entropy,sum_enthalpy,TS, sum_ZPE)

