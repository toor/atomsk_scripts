from snakemake.utils import min_version
import numpy as np
import datetime

min_version("7.0.0")

include: "jams.smk"

temperatures = [f'{x:.1f}' for x in range(10, 410, 10)]

cells = {
    "sc_100",
    "bcc_100",
    "bcc_110",
    "fcc_100",
    "fcc_111"
}

layer_i = 2
layer_f = 5

layers = range(layer_i, layer_f + 1)

rule all:
    input:
        expand("{lattice}/{layer}/{T}K/jams_mag.tsv",
            lattice=cells, layer=layers, T=temperatures),
        expand("{lattice}/{layer}/mag_vs_temp.png", lattice=cells, layer=layers)
