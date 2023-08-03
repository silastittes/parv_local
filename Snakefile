import re

my_scratch = "/scratch/stittes/"

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


REF, REGIONS = glob_wildcards("data/vcf/{ref}/vcffiles_split_raw/{regions}.vcf")

REGIONS = natural_sort(REGIONS)


#######
# mop #
#######

REF, SSP, POP = glob_wildcards("data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt")
mx_POP = [f"{i}--{j}" for i, j in  zip(SSP, POP)]
print(mx_POP)
mix_POP = list(set(list(filter(lambda a: a not in ['LR--random', 'Teo--random'], mx_POP))))
#print('\n'.join(mix_POP))

mSSP, mPOP, mCHROM, mSTART, mEND = glob_wildcards("data/mop/v5--{ssp}--{pop}--{chrom}--{start}--{end}.txt")
mop_files = expand("data/mop/v5--{ssp}--{pop}--{chrom}--{start}--{end}.bed", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)

chroms = list(set(mCHROM))

##########
##  PI  ##
##########

#SWITCH
pi_files = expand("data/angsd_pi/v5--{ssp}--{pop}--{chrom}--{start}--{end}.WINDOWBP_theta.thetasWindow.gz.pestPG", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
pi_files = [p.replace("WINDOW", w) for p in pi_files for w in ["1000", "100000", "1000000"]]

pi_full =  expand("data/angsd_pi/v5--{ssp}--{pop}.WINDOWBP_theta.txt", zip, pop = POP, ssp = SSP)
pi_full = [p.replace("WINDOW", w) for p in pi_full for w in ["1000", "100000", "1000000"]]

##########
##  VCF ##
##########

vcf_files = expand("data/angsd_vcf/v5--{ssp}--{pop}--{chrom}--{start}--{end}.vcf.mop.gz",  zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)


###############
## ngsRelate ##
###############

relate_files = expand("data/ngsRelate/v5--{ssp}--{pop}.ngsRelate.txt", zip, ssp = SSP, pop = POP)

##########
### MK ###
##########

#0: A,T -> G,C
#0: G,C -> A,T
#0: A,T -> A,T + G,C -> G,C

fold = ["0", "4"]
#WS, SW, SW+WS, ALL
nucs = ["AG,AC,TG,TC", "GA,GT,CA,CT", "AT,TA,GC,CG", "AT,AG,AC,TA,TG,TC,GA,GT,GC,CA,CT,CG"]

mk_sites = expand("data/angsd_sfs/v5--{ssp}--{pop}_foldFOLD_NUCS_sfs.txt",  zip, ssp = SSP, pop = POP)

mk_files = [m.replace("FOLD", f).replace("NUCS", n) for m in mk_sites for f in fold for n in nucs]

K_dict = {"Teo_k": [6], "LR_k": [5]}

###########
## RAiSD ##
###########

#
#dropping full pops for remaining analyses

raisd_corrected = expand("data/raisd/RAiSD_Report.v5--{ssp}--{pop}--{chrom}--{start}--{end}.corrected", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
raisd_outliers = expand("data/raisd/RAiSD_Report.v5--{ssp}--{pop}--{chrom}--{start}--{end}.corrected_block_outliers", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
raisd_merged = expand("data/raisd/v5--{ssp}--{pop}.corrected_block_outliers_merged.txt", zip, ssp = SSP, pop = POP)
print('\n'.join(raisd_outliers))

###########
## MUSHI ##
##########

mushi_out = expand("data/mushi/RAiSD_Report.v5--{ssp}--{pop}--msprime", zip, ssp = SSP, pop = POP)

#drop full PC files
raisd_corrected = [x for x in raisd_corrected if "LR--Palmar_Chico" not in x if "Teo--Palmar_Chico" not in x]
raisd_outliers = [x for x in raisd_outliers if "LR--Palmar_Chico" not in x if "Teo--Palmar_Chico" not in x]
raisd_merged = [x for x in raisd_merged if "LR--Palmar_Chico" not in x if "Teo--Palmar_Chico" not in x]
mushi_out = list(set([x for x in mushi_out if "LR--Palmar_Chico" not in x if "Teo--Palmar_Chico" not in x]))
pop_sweep_regions = expand("data/sweep_regions/sweeps.v5--{ssp}--{pop}--{chrom}--{start}--{end}.csv", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
pop_sweep_regions = list(set([x for x in pop_sweep_regions if "LR--Palmar_Chico" not in x if "Teo--Palmar_Chico" not in x]))

##########
## RDMC ##
##########

FREQ_POPS = [
    "v5--LR--Amatlan_de_Canas",
    "v5--LR--Crucero_Lagunitas",
    "v5--LR--Los_Guajes",
    "v5--LR--random1_Palmar_Chico",
    "v5--LR--San_Lorenzo",
    "v5--Teo--Amatlan_de_Canas",
    "v5--Teo--Crucero_Lagunitas",
    "v5--Teo--El_Rodeo",
    "v5--Teo--Los_Guajes",
    "v5--Teo--random1_Palmar_Chico",
    "v5--Teo--San_Lorenzo"
]


pop_freqs = [f"data/rdmc/freq/{popz}--{chrom}--{start}--{end}_freq.bed.gz" for chrom, start, end in list(set(zip(mCHROM, mSTART, mEND))) for popz in FREQ_POPS]

allpop_freqs = [f"data/rdmc/freq/v5--allpops--{chrom}--{start}--{end}_freq.txt.gz" for chrom, start, end in list(set(zip(mCHROM, mSTART, mEND)))]

#print('\n'.join(allpop_freqs))


############
## IBDseq ##
############

ibd_files = expand("data/ibdseq/v5--{ssp}--{pop}--{chrom}--{start}--{end}_r2maxR2.hbd", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
ibd_files = [f.replace("R2", r) for f in ibd_files for r in ["0.7", "0.4"]]

demography_ibd = expand("data/mushi/postsims/v5--{ssp}--{pop}_simREP.hbd", zip, ssp = SSP, pop = POP)
demography_ibd = [d.replace("REP", str(r)) for d in demography_ibd for r in range(10)]
#demography_ibd = [x for x in demography_ibd if "LR--Palmar_Chico" not in x if "Teo--Palmar_Chico" not in x]

all_files = []
for ref in ["v5"]:
  
    all_files.append(f"data/diplo/diplo_{ref}.fa.gz")
    all_files.append(f"data/trip/trip_{ref}.fa.gz")
    all_files.append(f"data/lux/{ref}--lux--lux_nuctable.txt")
    all_files.append(f"data/diplo/{ref}--diplo--diplo_nuctable.txt")
    all_files.append(f"data/diplo/{ref}--diplo--diplo.mafs.gz")
    all_files.append(f"data/anc/{ref}_anc.fa")
    all_files.append(f"data/refs/{ref}/{ref}_FOLD")
    all_files.append(f"data/anc/{ref}_anc_FOLD")
    all_files.append(f"data/angsd_treemix/{ref}_treemix.treeout.gz")

    for ssp in ["Teo", "LR"]:
        ref_dict = K_dict[f'{ssp}_k']
        all_files.append(f"data/beagle/{ref}--{ssp}.beagle.gz")
        for k in ref_dict:
            all_files.append(f"data/ngsAdmix/{ref}--{ssp}_K{k}.qopt")
            all_files.append(f"data/ngsAdmix/v5_{ssp}_{k}_thin1M_random10_PalmarChico.qopt")
           
        all_files.append(f"data/refs/{ref}/{ref}_FOLD")
        all_files.append(f"data/trip/trip_{ref}_FOLD")
        all_files.append(f"data/beagle/v5--{ssp}_100thin_random10_PalmarChico.beagle.gz")
        all_files.append(f"data/beagle/v5--{ssp}_100thin_PalmarChicoONLY.beagle.gz")
        for i in range(2,6,1):
            all_files.append(f"data/ngsAdmix/v5_{ssp}_{i}_thin1M_PalmarChicoONLY.qopt")    
 
rule all:
    input:
        all_files, 
        "data/diplo/v5--diplo--diplo.mafs.gz",
        mop_files,
        pi_files,
        pi_full,
        mk_files,
        relate_files,
        vcf_files,
        ibd_files,
        demography_ibd,
        mushi_out,
        raisd_corrected,
        #raisd_outliers,
        #raisd_merged,
        #"data/raisd/v5--allpops--shared_outliers.txt",
        pop_sweep_regions,
        "data/sweep_regions/v5--allpops--shared_outliers.txt",
        pop_freqs,
        allpop_freqs

include: "rules/process_raw.smk"
include: "rules/mop.smk"
include: "rules/outgroup.smk"
include: "rules/angsd.smk"
include: "rules/mk.smk"
include: "rules/ngsRelate.smk"
include: "rules/treemix.smk"
include: "rules/angsd_vcf.smk"
include: "rules/demography.smk"
include: "rules/rdmc.smk"


