import numpy as np
import matplotlib.pyplot as plt

sep = '    '

input_file = snakemake.input
output_file = snakemake.output

data = np.genfromtxt(input_file, sep)

sorted_data = data[data[:,0].argsort()]

plt.figure()

plt.title("blah")
plt.xlabel("Temperature")
plt.ylabel("Magnetisation")

temps = data[:,0]
mags = data[:,1]

plt.plot(temps, mags)
plt.savefig(output_file)
