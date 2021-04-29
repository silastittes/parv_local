#input sfs, corrected raisd file, file prefix  

#fit demographic model
#print obs and exp sfs to figure
#print prefix, parameters, etc. to data file for alter concactenation
#run msprime under the demographic model


import mushi
import msprime
import tskit
import argparse
import pandas as pd
import numpy as np
import os
from itertools import tee


parser = argparse.ArgumentParser(
    prog = "Simulate demography using mushi, generate data using msprime, run ibdSeq on simulations",    
    description=""
)

parser.add_argument('-s', '--sfs', type=str, nargs='+', required = True,
            help='One or more input files containing the unfolded sfs on the first line, white space delimited.')

parser.add_argument('-p', '--prefix', type=str, required = True,
            help='String to add to output files. If a path is included, raisd files will be named my the final field after splitting on back slashes.')


args = parser.parse_args()

sfs_list = []

for sfs_file in args.sfs:
    with open(sfs_file) as f:
        sfs = f.readlines()[0].split()
        sfs_chr = np.array([float(s) for s in sfs])
        sfs_list.append(sfs_chr)

sfs_array = np.asarray(sfs_list)
sfs_all = sfs_array.sum(axis=0)


print(sfs_all)


def sim_hbd(sfs_all, prefix, mu = 3e-8, c = 1.6e-8, bps = 5e6, reps = 10, strong_reg = True):
    
    t = np.logspace(np.log10(1), np.log10(1000000), 50)
    ksfs = mushi.kSFS(np.array(sfs_all[1:-1]))
    
    if strong_reg:
        ksfs.infer_history(t, mu0 = mu*sum(sfs_all), infer_mu=False, folded = False,
                           alpha_tv=1e4, alpha_spline=1e4, alpha_ridge = 1e-1,
                           tol=1e-12, verbose=True, max_iter = 5000)
    else:
        ksfs.infer_history(t, mu0 = mu*sum(sfs_all), infer_mu=False, folded = False,
                           tol=1e-12, verbose=True, max_iter = 5000)

    Nt = ksfs.eta.vals
    T = ksfs.eta.change_points

    N_0 = Nt[0] / 2
    msp_N = Nt[1:] / 2
    msp_T = T
    nsamp = len(sfs_all)-1

    population_configurations = [msprime.PopulationConfiguration(sample_size = nsamp, initial_size = N_0)]
    demography_list = [msprime.PopulationParametersChange(time = t, initial_size = n) for t, n, in zip(msp_T, msp_N)]

    ts_mushi = msprime.simulate(
        population_configurations=population_configurations, 
        demographic_events = demography_list, 
        length = bps, 
        recombination_rate = c, 
        mutation_rate = mu, 
        num_replicates = reps
    )
    
    
    for i,ts in enumerate(ts_mushi):
        ts = msprime.mutate(ts, rate=mu, model=msprime.InfiniteSites(alphabet=msprime.NUCLEOTIDES))
        with open(f"{prefix}_sim{i}.vcf", "w") as vcf_file:
            ts.write_vcf(vcf_file, ploidy=2)
        
        os.system(f"java -jar src/ibdseq.r1206.jar gt={prefix}_sim{i}.vcf out={prefix}_sim{i} nthreads=2 r2max=0.7")
        

sim_hbd(sfs_all, f"{args.prefix}", mu = 3e-8, c = 1.6e-8, bps = 1e7, reps = 10)

