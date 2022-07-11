import numpy as np
import subprocess

jams_cfg = snakemake.input[0]
jams_h5 = snakemake.input[1]

cmd = []

cmd.append(exe)
cmd.append(f"--output=\"{output}\"")
cmd.append("--name")
cmd.append(name)
cmd.append(f" \"lattice : {{spins=\\\"{jams_h5}\\\";}};\" ")
cmd.append(jams_cfg)

subprocess.run(cmd)
