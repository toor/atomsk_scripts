from snakemake.utils import min_version
import numpy as np
import datetime

min_version("7.0.0")

include: "jams.smk"

temperatures = [f'{x:.1f}' for x in range(10, 510, 10)]
#moment_values = [f'{x:.1f}' for x in [3.0, 5.0, 7.0]]
thermostat = 'langevin-bose-gpu'
#repeats = [f'{x:.1f}' for x in np.arange(1,21,1)]

cells = {
    "sc_100",
    "bcc_100",
    "bcc_110",
    "fcc_100",
    "fcc_111"
}

layer_i = 1
layer_f = 5

layers = range(layer_i, layer_f + 1)

rule all:
    input:
        expand("{lattice}/{layer}/{T}K/jams_mag.tsv",
            lattice=cells, layer=layers, T=temperatures)
