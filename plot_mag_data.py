import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, m_e, hbar, e

bohr = e*hbar/(2*m_e)

sep = '\t'

input_file = str(snakemake.input)
output_file = str(snakemake.output)

lattice = snakemake.params[0]
layer = int(snakemake.params[1])
constant = float(snakemake.params[2])

print(layer)

d_planes = {
    "sc_100": constant,
    "bcc_100": constant/2,
    "bcc_110": constant/np.sqrt(2),
    "fcc_100": constant/2,
    "fcc_111": constant/np.sqrt(3)
}

# Bloch's law in thin films
def bloch(T):
    # volume of Wigner-Seitz unit cell (1st BZ)
    v_ws = 2*np.pi/np.power(constant, 3)
    L_z = layer*d_planes[lattice]
    
    # TODO make all of these parameters in the Snakefile, to also pass to e.g. jinja
    # Exchange
    J = 2e-21
    # TODO: Which value to use for the total spin? How is this calculated? am I just being thick?
    S = 0.5
    # Spinwave stiffness
    D = 2*J*S*(constant**2)

    prefactor = v_ws*k*T/(4*np.pi*S*L_z*D) 

    array = np.zeros((layer, 1))
    for i in range(0, layer):
        k_z = i*np.pi/(constant*(layer - 1))

        exp = np.exp((-D/(k*T))*np.power(k_z, 2))
        print("exp = " + str(exp))
        array[i] = np.log(1 - exp)
    
    m_z = 1 + prefactor*np.sum(array)

    return m_z

data = []

with open(input_file) as f:
    for line in f:
        l = [float(i) for i in line.strip().split(sep)]
        data.append(l)
data = np.array(data)

plt.figure()

title = "|M_z| vs. Temperature for "

plt.title("|M_z| vs. Temperature.")
plt.xlabel("Temperature (K)")
plt.ylabel("Magnetisation (T)")

temps = data[:,0]
mags = data[:,1]

bloch_data = np.array([bloch(T) for T in temps])

plt.plot(temps, mags, label='JAMS data')
plt.plot(temps, bloch_data, label='Bloch law')
plt.savefig(output_file)
