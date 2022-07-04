import numpy as np
import matplotlib.pyplot as plt

sep = '\t'

input_file = snakemake.input[0]
output_file = snakemake.output[0]


#data = np.genfromtxt(input_file, sep)

data = []
with open(input_file) as f:
    for line in f:
        l = [float(i) for i in line.strip().split(sep)]
        data.append(l)

data = np.array(data)

plt.figure()

plt.title("blah")
plt.xlabel("Temperature")
plt.ylabel("Magnetisation")

# Extract 1st and 2nd columns
temps = data[:,0]
mags = data[:,1]

plt.plot(temps, mags)
plt.savefig(output_file)
