#input sfs, corrected raisd file, file prefix  

#fit demographic model
#print obs and exp sfs to figure
#print prefix, parameters, etc. to data file for alter concactenation
#run msprime under the demographic model


import mushi
import msprime
import tskit
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from itertools import tee


parser = argparse.ArgumentParser(
    prog = "Simulate demography using mushi, generate data using msprime to compare simulated and observed summary statistics.",    
    description=""
)

parser.add_argument('-i', '--ID', type=str, required = True,
            help='population ID to include in output files where appropriate.')

parser.add_argument('-s', '--sfs', type=str, nargs='+', required = True,
            help='One or more input files containing the unfolded sfs on the first line, white space delimited.')

parser.add_argument('-p', '--prefix', type=str, required = True,
            help='String to add to output files. If a path is included, raisd files will be named my the final field after splitting on back slashes.')

parser.add_argument('-n', '--nsims', type=int, required = True,
            help='The number of simulations do generate with msprime under the mushi demography.')

parser.add_argument('-b', '--bps', type=int, required = True,
            help='The number of basepairs to simulate with msprime.')


args = parser.parse_args()



mu = 1e-7
c = 1.6e-8

first = True
for sfs_file in args.sfs:
    with open(sfs_file) as f:
        sfs = f.readlines()[0].split()
        sfs = np.array([float(s) for s in sfs])
    if first:
        sfs = sfs
        first = False
    else:
        sfs += sfs

print(sfs)

t = np.logspace(np.log10(1), np.log10(500000), 100)
ksfs = mushi.kSFS(np.array(sfs[1:-1]))
ksfs.infer_history(t, mu0 = mu*sum(sfs), infer_mu=False, folded = False, verbose=True, max_iter = 5000,  tol=1e-11)

ksfs.plot_total()
plt.xscale('log')
plt.yscale('log')
plt.savefig(f"{args.prefix}_sfs.png")

ksfs.eta.plot()
plt.savefig(f"{args.prefix}_demography.png")


Nt = ksfs.eta.vals
T = ksfs.eta.change_points

nsims = args.nsims
bps = args.bps
ns = len(sfs)-1
N_0 = Nt[0] / 2
msp_N = Nt[1:] / 2
msp_T = T  


population_configurations = [msprime.PopulationConfiguration(sample_size = ns, initial_size = N_0)]
demography_list = [msprime.PopulationParametersChange(time = t, initial_size = n) for t, n, in zip(msp_T, msp_N)]

ts_mushi = msprime.simulate(
    population_configurations=population_configurations, 
    demographic_events = demography_list, 
    #Ne = N_0, sample_size = ns,
    length = bps, 
    recombination_rate = c, 
    mutation_rate = mu, 
    num_replicates = nsims
)

ts_mushi, ts_mushi_copy = tee(ts_mushi)

pop = args.ID
with open(f"{args.prefix}_demography.txt", "w") as demography_file:
    for t, n in zip(msp_T, msp_N):
        print(f"{t}\t{n}\t{pop}", file = demography_file)
 
with open(f"{args.prefix}_stats.txt", "w") as stat_file:
    for ts in ts_mushi:
        print(f"{ts.diversity().item()}\t{ts.Tajimas_D().item()}\t{pop}", file = stat_file)


with open(f"{args.prefix}_ms.txt", "w") as ms_file:
    tskit.write_ms(ts_mushi_copy, output = ms_file)

