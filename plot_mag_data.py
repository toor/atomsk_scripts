import numpy as np
import matplotlib.pyplot as plt

sep = '\t'

input_file = str(snakemake.input)
output_file = str(snakemake.output)

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
bloch_mags = data[:,2]

plt.plot(temps, mags, label='JAMS data')
plt.plot(temps, bloch_mags, label='Bloch law')
plt.legend()

plt.savefig(output_file)
