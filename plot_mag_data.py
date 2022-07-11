import numpy as np
import matplotlib.pyplot as plt

sep = '\t'

mag_input_file = str(snakemake.input[0])
sus_input_file = str(snakemake.input[1])

mag_output_file = str(snakemake.output[0])
sus_output_file = str(snakemake.output[1])

mag_data = np.genfromtxt(mag_input_file, skip_header=1)
sus_data = np.genfromtxt(sus_input_file, skip_header=1, usecols=1)

temps = mag_data[:,0]
m_z = mag_data[:,1]
m = mag_data[:,2]
m_bloch = mag_data[:,3]

plt.figure()

title = "|M_z| vs. Temperature for "

plt.title("m_z vs. Temperature.")
plt.xlabel("Temperature (K)")
plt.ylabel("Magnetisation (T)")

plt.plot(temps, m_z, label='JAMS M_z')
plt.plot(temps, m, label='JAMS |M|')
plt.plot(temps, m_bloch, label='Bloch law')
plt.legend()

plt.savefig(mag_output_file)

plt.figure()

plt.title("Susceptibility vs. temperature")
plt.xlabel("Temperature (K)")
plt.ylabel("Temperature-normalised susceptibility (T^2)")

plt.plot(temps, sus_data)

plt.savefig(sus_output_file)
