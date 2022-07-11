import os
import numpy as np
import scipy
import pandas as pd
import matplotlib.pyplot as plt
import re
from pymbar import timeseries
from scipy.constants import k, m_e, hbar, e

TEMP_REGEX = r'[+-]?([0-9]*[.])[0-9]+'

# use 4 space delimiter

sep = '\t'

in_files = snakemake.input
mag_out_file = snakemake.output[0]
sus_out_file = snakemake.output[1]


lattice = snakemake.params[0]
layer = int(snakemake.params[1])
constant = float(snakemake.params[2])
app_field = float(snakemake.params[3])
repeat = int(snakemake.params[4])

# number of atoms in the system
N = repeat*repeat*layer
l_x = (repeat - 1)*constant
l_y = l_x
l_z = (layer - 1)*constant

V = l_x*l_y*l_z
v_ws = V/N

d_planes = {
    "sc_100": constant,
    "bcc_100": constant/2,
    "bcc_110": constant/np.sqrt(2),
    "fcc_100": constant/2,
    "fcc_111": constant/np.sqrt(3)
}

# Bloch's law in thin films
def bloch(T):
    bohr = (e*hbar)/(2*m_e)
    g = 2
    
    # volume of Wigner-Seitz unit cell (1st BZ)
    #v_ws = 2*np.pi/np.power(constant, 3)
    L_z = layer*d_planes[lattice]
    
    # TODO make all of these parameters in the Snakefile, to also pass to e.g. jinja
    # Exchange
    J = 2e-21
    S = 1.25
    # Spinwave stiffness
    D = 2*J*S*(constant**2)

    prefactor = v_ws*k*T/(4*np.pi*S*L_z*D) 

    array = np.zeros((layer, 1))
    for i in range(0, layer):
        k_z = i*np.pi/(constant*(layer - 1))

        exp = np.exp(-1/(k*T)*(g*bohr*app_field + D*np.power(k_z, 2)))
        print("exp = " + str(exp))
        array[i] = np.log(1 - exp)
    
    m_z = 1 + prefactor*np.sum(array)

    return m_z

with open(mag_out_file) as f:
    w_str = "T" + sep + "M_z" + sep + "M" + sep + "M_bloch" + "\n"
    f.write(w_str)
    f.close()
with open(sus_out_file) as f:
    w_str = "T" + sep + "Sus" + "\n"

for in_file in in_files:
    # Use all of the magnetisation data, since we start from a thermalised state
    m_z_data = np.genfromtxt(in_file, skip_header=1, usecols=7)
    m_data = np.genfromtxt(in_file, skip_header=1, usecols=8)

    m_z = np.mean(m_z_data)
    mag = np.mean(m_data)
    sus = np.var(m_data) - mag**2

    m = re.search(TEMP_REGEX, in_file)
    temp = float(in_file[m.start():m.end()])
    mag_bloch = bloch(temp)
    diff = np.abs(mag - mag_bloch)

    mag_str = str(temp) + sep + str(round(m_z, 3)) + sep + str(round(mag, 3)) + sep + str(round(mag_bloch, 3)) + "\n"
    sus_str = str(temp) + sep + str(round(sus, 3)) + "\n"

    with open(str(mag_out_file)) as f:
        f.write(mag_str)
        f.close()
    with open(str(sus_out_file)) as f:
        f.write(sus_str)
        f.close()

