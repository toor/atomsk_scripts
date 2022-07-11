from snakemake.utils import min_version
import numpy as np
import datetime

min_version("7.0.0")

include: "jams.smk"

temperatures = [f'{x:.1f}' for x in range(25, 550, 25)]
#temperatures = [50.0, 100.0, 200.0, 300.0, 400.0]

cells = {
    "sc_100"
    #"bcc_100",
    #"bcc_110",
    #fcc_100"
    #"fcc_111"
}

layer_i = 2
layer_f = 5


layers = range(layer_i, layer_f + 1)
#layers = 3
#temp = 50.0

rule all:
    input:
        expand("{lattice}/{layer}/{T}K/jams_mag.tsv",
            lattice=cells, layer=layers, T=temperatures),
        expand("{lattice}/{layer}/mag_vs_temp.png", lattice=cells, layer=layers),
        expand("{lattice}/{layer}/sus_vs_temp.png", lattice=cells, layer=layers),
        expand("figures/magnon_spectrum_{lattice}_{layer}_{T}K.pdf",
            lattice=cells, layer=layers, T=temperatures)
