#input sfs, corrected raisd file, file prefix  

#fit demographic model
#print obs and exp sfs to figure
#print prefix, parameters, etc. to data file for alter concactenation
#run msprime under the demographic model


import mushi
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

parser = argparse.ArgumentParser(
    prog = "Simulate demography using mushi, generate data using msprime to compare simulated and observed summary statistics.",    
    description=""
)

parser.add_argument('-i', '--ID', type=str, required = True,
            help='population ID to include in output files where appropriate.')

parser.add_argument('-s', '--sfs', type=str, required = True,
            help='Input file containing the unfolded sfs on the first line, white space delimited.')

parser.add_argument('-p', '--prefix', type=str, required = True,
            help='String to add to output files. If a path is included, raisd files will be named my the final field after splitting on back slashes.')

parser.add_argument('-n', '--nsims', type=int, required = True,
            help='The number of simulations do generate with msprime under the mushi demography.')

parser.add_argument('-b', '--bps', type=int, required = True,
            help='The number of basepairs to simulate with msprime.')


args = parser.parse_args()



mu = 1e-7
with open(args.sfs) as f:
    sfs = f.readlines()[0].split()
    sfs = [float(s) for s in sfs]

t = np.logspace(np.log10(1), np.log10(1000000), 100)
ksfs = mushi.kSFS(np.array(sfs[1:-1]))
ksfs.infer_history(t, mu0 = mu*sum(sfs), infer_mu=False, folded = False, verbose=True, max_iter = 10000)

ksfs.plot_total()
plt.xscale('log')
plt.yscale('log')
plt.savefig(f"{args.prefix}_sfs.png")

ksfs.eta.plot()
plt.savefig(f"{args.prefix}_demography.png")


Nt = ksfs.eta.vals
T = ksfs.eta.change_points

N_0 = Nt[0]
ms_N = Nt[1:]/N_0
ms_T = T/(4*N_0)

pop = args.ID

demography_file = open(f"{args.prefix}_demography.txt", "w")
#print(f"t\tn\tpop")
for t, n in zip(ms_T, ms_N):
    print(f"{t}\t{n}\t{pop}", file = demography_file)
 
en_string = ' '.join([f"--size-change {t} {n}" for t, n in zip(ms_T, ms_N)])

ns = len(sfs)-1
bps = args.bps
nsims = args.nsims
c = 1.6e-8
theta = 4 * N_0 * 3e-8 * bps
rho = 4 * N_0 * c * bps

ms_string = f"mspms {ns} {nsims} --mutation-rate {theta} --recombination {rho} {bps} {en_string} > {args.prefix}_msprime.txt"
os.system(ms_string)
msp_command_file = open(f"{args.prefix}_mspms_command.txt", "w")
print(ms_string, file = msp_command_file)

