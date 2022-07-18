import os
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import re
from pymbar import timeseries as ts
from scipy.constants import k, m_e, hbar, e
from scipy.special import zeta, gamma


in_files = snakemake.input
mag_out_file = snakemake.output[0]
sus_out_file = snakemake.output[1]

lattice = snakemake.params[0]
layer = int(snakemake.params[1])
constant = float(snakemake.params[2])
app_field = float(snakemake.params[3])
repeat = int(snakemake.params[4])

d_planes = {
    "sc_100": constant,
    "bcc_100": constant/2,
    "bcc_110": constant/np.sqrt(2),
    "fcc_100": constant/2,
    "fcc_111": constant/np.sqrt(3)
}

TEMP_REGEX = r'[+-]?([0-9]*[.])[0-9]+'

# use 4 space delimiter

sep = '\t'
# number of atoms in the system
N = repeat*repeat*layer
l_x = (repeat - 1)*constant
l_y = l_x
l_z = (layer - 1)*d_planes[lattice]

V = l_x*l_y*l_z
v_ws = V/N

# Bloch's law in thin films
def bloch_thinfilms(T):
    # bohr magneton, mu_B = 5.7883818012×10−5 eV*T^-1
    bohr = (e*hbar)/(2*m_e)
    g = 2
 
    # TODO make all of these parameters in the Snakefile, to also pass to e.g. jinja
    # Exchange
    J = 2e-21
    S = int(snakemake.params.get("mu_s"))/2
    # Spinwave stiffness
    D = 2*J*S*(constant**2)

    prefactor = v_ws*k*T/(4*np.pi*S*l_z*D) 

    array = np.zeros((layer, 1))
    for i in range(0, layer):
        k_z = i*np.pi/(constant*(layer - 1))

        exp = np.exp(-1/(k*T)*(g*bohr*app_field + D*np.power(k_z, 2)))
        array[i] = np.log(1 - exp)
    
    m_z = 1 + prefactor*np.sum(array)

    return m_z

def bloch_bulk(T):
    bohr = (e*hbar)/(2*m_e)
    g = 2

    S = int(snakemake.params.get("mu_s"))/2
    J = 2e-21
    D = 2*J*S*(constant**2) 

    # saturation magnetisation at T=0K
    M_zs = g*bohr*S/v_ws
    prefactor = g*bohr*zeta(3/2)*gamma(3/2)/(M_zs*(4*(np.pi**2)))

    m_z = 1 - prefactor*np.power((k*T)/D, 3/2)

    return m_z

with open(mag_out_file, 'a') as f:
    w_str = "T" + sep + "M_z" + sep + "M" + sep + "M_bloch (thin film)" + sep + "M_bloch (bulk)" + "\n"
    f.write(w_str)
    f.close()
with open(sus_out_file, 'a') as f:
    w_str = "T" + sep + "Sus" + "\n"
    f.write(w_str)
    f.close()

for in_file in in_files:
    # Use all of the magnetisation data, since we start from a thermalised state
    m_z_data = np.genfromtxt(in_file, skip_header=1, usecols=7)
    m_data = np.genfromtxt(in_file, skip_header=1, usecols=8)
    
    m = re.search(TEMP_REGEX, in_file)
    temp = float(in_file[m.start():m.end()])
    mag_bloch_thin = bloch_thinfilms(temp)
    mag_bloch_bulk = bloch_bulk(temp)
    
    m_z = np.mean(m_z_data[4000:])
    mag = np.mean(m_data[4000:])
    sus = np.var(m_data[4000:])

    t0_mz = ts.detectEquilibration(m_z_data)[0]
    t0_m = ts.detectEquilibration(m_data)[0]
    
    #print("Using an MBAR estimator, detected the end of the equilibration period for m_z at timestep index " + str(t0_mz) + " at T=" + str(temp))
    print("Using an MBAR estimator, detected the end of the equilibration period for m at timestep index " + str(t0_m) + " at T=" + str(temp))

    mag_str = str(temp) + sep + str(round(m_z, 3)) + sep + str(round(mag, 3)) + sep + str(round(mag_bloch_thin, 3)) + sep + str(round(mag_bloch_bulk, 3)) + "\n"
    sus_str = str(temp) + sep + str(round(sus, 3)) + "\n"

    with open(str(mag_out_file), 'a') as f:
        f.write(mag_str)
        f.close()
    with open(str(sus_out_file), 'a') as f:
        f.write(sus_str)
        f.close()
