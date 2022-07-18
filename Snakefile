from snakemake.utils import min_version
import numpy as np
import datetime

min_version("7.0.0")

layers = range(3, 7)

include: "jams.smk"

temperatures = [f'{x:.1f}' for x in range(10, 710, 10)]
magnon_temps = [f'{x:.1f}' for x in [10, 25, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700]]

cells = {
    #"sc_100"
    "bcc_100",
    "bcc_110",
    "fcc_100",
    "fcc_111"
}

rule all:
    input:
        expand("results/{lattice}/{layer}/{T}K/jams_mag.tsv",
            lattice=cells, layer=layers, T=temperatures),
        expand("results/{lattice}/{layer}/mag_vs_temp.png", lattice=cells, layer=layers),
        expand("results/{lattice}/{layer}/sus_vs_temp.png", lattice=cells, layer=layers),
        expand("figures/{lattice}/{layer}/magnon_{T2}K/spectrum.pdf",
            lattice=cells, layer=layers, T2=magnon_temps)
