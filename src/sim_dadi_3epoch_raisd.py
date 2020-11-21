#input sfs, corrected raisd file, file prefix  

#fit demographic model
#print obs and exp sfs to figure
#print prefix, parameters, etc. to data file for alter concactenation
#run msprime under the demographic model

import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    prog = "tidy_vcf",    
    description="Given a vcf file, produces a tidy versions of sites and genotypes data, in efforts to make it easier to calculate summary statitics and visualizations."
)


parser.add_argument('-s', '--sfs', type=str, required = True,
            help='Input file containing the unfolded sfs on the first line, white space delimited.')

parser.add_argument('-p', '--prefix', nargs="?", type=str, required = True,
            help='String to add to output files. If a path is included, raisd files will be named my the final field after splitting on back slashes.')

parser.add_argument('-n', '--nsims', nargs="?", type=int, required = True,
            help='The number of simulations do generate with msprime under the ML 2-epoch demography.')

parser.add_argument('-b', '--bps', nargs="?", type=int, required = True,
            help='The number of basepairs to simulate with msprime.')


args = parser.parse_args()


def three_epoch(params, ns, pts):
    """
    params = (nuB,nuF,TB,TF)
    ns = (n1,)

    nuB: Ratio of bottleneck population size to ancient pop size
    nuF: Ratio of contemporary to ancient pop size
    TB: Length of bottleneck (in units of 2*Na generations) 
    TF: Time since bottleneck recovery (in units of 2*Na generations) 

    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,TB,TF = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, TB, nuB)
    phi = Integration.one_pop(phi, xx, TF, nuF)

    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs



with open(args.sfs) as f:
    sfs = f.readlines()[0].split()
    sfs = [float(s) for s in sfs]

fs = dadi.Spectrum(np.round(sfs, 0))

params = (1, 1, 1, 1)
ns = fs.sample_sizes
pts = [200]

my_extrap_func = Numerics.make_extrap_log_func(three_epoch)

model = my_extrap_func(params, ns, pts)

ll_model = dadi.Inference.ll_multinom(model, fs)

theta = dadi.Inference.optimal_sfs_scaling(model, fs)

lower_bound, upper_bound = [1e-2, 1e-2, 1e-2, 1e-2], [1e2, 1e2, 1e2, 1e2]
p0 = dadi.Misc.perturb_params(params, lower_bound=lower_bound, upper_bound=upper_bound)
dadi_opt = dadi.Inference.optimize_log(p0, fs, my_extrap_func, pts, verbose=1, maxiter=10000, lower_bound=lower_bound, upper_bound=upper_bound)

#p0 = dadi.Misc.perturb_params(params)
#dadi_opt = dadi.Inference.optimize_log(p0, fs, my_extrap_func, pts, verbose=1, maxiter=1000)


plt.plot(sfs[1:-1]/theta)
plt.plot(my_extrap_func(params, ns, pts)[1:-1])
plt.xlabel("Allele frequency") 
plt.ylabel("Frequency")
plt.title(f"{args.prefix}")
plt.savefig(f"{args.prefix}_sfs.png")

ns = ns
bps = args.bps
mu = 3e-8
c = 1.6e-8

nu_B, nu_F, T_B, T_F = dadi_opt
theta_sc = theta / sum(sfs)
N_anc = theta_sc/(4*mu)
N_0 = N_anc * nu_F
N_B = N_anc * nu_B


stat_file =  open(f"{args.prefix}_dadi_stat.txt", "w")

print("pop\tnu_F\tnu_B\tT_B\tT_F\ttheta\tN_0\tN_A\tN_B", file = stat_file)
print(f"{args.prefix}\t{nu_F}\t{nu_B}\t{T_B}\t{T_F}\t{theta_sc}\t{N_0}\t{N_anc}\t{N_B}", file = stat_file)

msp_command = f"mspms {ns[0]} {args.nsims} --mutation-rate {4 * N_0 * mu * bps} --recombination {4 * N_0 * c * bps} {bps} --size-change {T_F/2} {N_B/N_0} --size-change {(T_F+T_B)/2} {N_anc/N_0} > {args.prefix}_msprime.txt"
os.system(msp_command)

msp_command_file = open(f"{args.prefix}_mspms_command.txt", "w")
print(msp_command, file = msp_command_file)
