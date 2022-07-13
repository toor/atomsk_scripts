from snakemake.utils import min_version
import numpy as np
import datetime

min_version("7.0.0")

layers = range(2, 7)

include: "jams.smk"

temperatures = [f'{x:.1f}' for x in range(50, 550, 50)]

cells = {
    #"sc_100"
    "bcc_100",
    "bcc_110"
    #fcc_100"
    #"fcc_111"
}

rule all:
    input:
        expand("{lattice}/{layer}/{T}K/jams_mag.tsv",
            lattice=cells, layer=layers, T=temperatures),
        expand("{lattice}/{layer}/mag_vs_temp.png", lattice=cells, layer=layers),
        expand("{lattice}/{layer}/sus_vs_temp.png", lattice=cells, layer=layers),
        #expand("figures/magnon_spectrum_{lattice}_{layer}_{T}K.pdf",
        #    lattice=cells, layer=layers, T=temperatures)
