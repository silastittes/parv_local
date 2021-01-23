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
#data/bamlist/til11--LR--Amatlan_de_Canas__bamlist.txt

REF, SSP, POP = glob_wildcards("data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt")


mx_POP = [f"{i}--{j}" for i, j in  zip(SSP, POP)]
print(mx_POP)
mix_POP = list(set(list(filter(lambda a: a not in ['LR--random', 'Teo--random'], mx_POP))))
#mix_POP.append("trip--trip")
print('\n'.join(mix_POP))


#SWITCH
#mREF, mSSP, mPOP, mCHROM, mSTART, mEND = glob_wildcards("data/mop/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.txt")
#mop_files = expand("data/mop/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.bed", zip, ref = mREF, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)

mSSP, mPOP, mCHROM, mSTART, mEND = glob_wildcards("data/mop/v5--{ssp}--{pop}--{chrom}--{start}--{end}.txt")
mop_files = expand("data/mop/v5--{ssp}--{pop}--{chrom}--{start}--{end}.bed", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)

chroms = list(set(mCHROM))

#mop_final = expand("data/mop/{ref}--{ssp}--{pop}_all.bed", zip, ref = mREF, ssp = mSSP, pop = mPOP)


##########
##  PI  ##
##########

#SWITCH
#pi_files = expand("data/angsd_pi/{ref}--{ssp}--{pop}.{win}kb_theta.thetasWindow.gz.pestPG", zip, ref = REF, pop = POP, ssp = SSP)
#pi_files = expand("data/angsd_pi/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.WINDOWBP_theta.thetasWindow.gz.pestPG", zip, ref = mREF, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
pi_files = expand("data/angsd_pi/v5--{ssp}--{pop}--{chrom}--{start}--{end}.WINDOWBP_theta.thetasWindow.gz.pestPG", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
pi_files = [p.replace("WINDOW", w) for p in pi_files for w in ["1000", "100000", "1000000"]]

pi_full =  expand("data/angsd_pi/v5--{ssp}--{pop}.WINDOWBP_theta.txt", zip, pop = POP, ssp = SSP)
pi_full = [p.replace("WINDOW", w) for p in pi_full for w in ["1000", "100000", "1000000"]]

##########
##  VCF ##
##########

#SWITCH
#data/angsd_vcf/{ref}--{ssp}--{pop}.vcf
#vcf_files = expand("data/angsd_vcf/{ref}--{ssp}--{pop}.vcf.gz", zip, ref = REF, pop = POP, ssp = SSP)
#vcf_files = expand("data/angsd_vcf/v5--{ssp}--{pop}.vcf.gz", zip, ref = REF, pop = POP, ssp = SSP)


###############
## ngsRelate ##
###############

#SWITCH
#relate_files = expand("data/ngsRelate/{ref}--{ssp}--{pop}.ngsRelate.txt", zip, ref = REF, ssp = SSP, pop = POP)
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
#nucs = ["AT,AG,AC,TA,TG,TC,GA,GT,GC,CA,CT,CG"]


#SWITCH
#mk_sites = expand("data/angsd_sfs/{ref}--{ssp}--{pop}_foldFOLD_NUCS_sfs.txt",  zip, ref = REF, ssp = SSP, pop = POP)
mk_sites = expand("data/angsd_sfs/v5--{ssp}--{pop}_foldFOLD_NUCS_sfs.txt",  zip, ssp = SSP, pop = POP)

mk_files = [m.replace("FOLD", f).replace("NUCS", n) for m in mk_sites for f in fold for n in nucs]

#print('\n'.join(mk_files))

K_dict = {"Teo_k": [6], "LR_k": [5]}

###########
## RAiSD ##
###########

#data/raisd/RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{region}.corrected
#SWITCH
#raisd_files = expand("data/raisd/RAiSD_Report.{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.bed", zip, ref = mREF, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)

#raisd_files = expand("data/raisd/RAiSD_Report.v5--{ssp}--{pop}--{chrom}--{start}--{end}.txt", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
raisd_corrected = expand("data/raisd/RAiSD_Report.v5--{ssp}--{pop}--{chrom}--{start}--{end}.corrected", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
raisd_outliers = expand("data/raisd/RAiSD_Report.v5--{ssp}--{pop}--{chrom}--{start}--{end}.corrected_block_outliers", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
raisd_merged = expand("data/raisd/v5--{ssp}--{pop}.corrected_block_outliers_merged.txt", zip, ssp = SSP, pop = POP)
print('\n'.join(raisd_outliers))



##########
## RDMC ##
##########

FREQ_POPS = [
    "v5--LR--Amatlan_de_Canas",
    "v5--LR--Crucero_Lagunitas",
    "v5--LR--Los_Guajes",
    "v5--LR--random1_Palmar_Chico",
    "v5--LR--random2_Palmar_Chico",
    "v5--LR--San_Lorenzo",
    "v5--Teo--Amatlan_de_Canas",
    "v5--Teo--Crucero_Lagunitas",
    "v5--Teo--El_Rodeo",
    "v5--Teo--Los_Guajes",
    "v5--Teo--random1_Palmar_Chico",
    "v5--Teo--random2_Palmar_Chico",
    "v5--Teo--San_Lorenzo"
]


pop_freqs = [f"data/rdmc/freq/{popz}--{chrom}--{start}--{end}_freq.bed.gz" for chrom, start, end in list(set(zip(mCHROM, mSTART, mEND))) for popz in FREQ_POPS]

allpop_freqs = [f"data/rdmc/freq/v5--allpops--{chrom}--{start}--{end}_freq.txt.gz" for chrom, start, end in list(set(zip(mCHROM, mSTART, mEND)))]

print('\n'.join(allpop_freqs))


##########
## DADI ##
##########

#SWTICH
#dadi_out = expand("data/dadi/RAiSD_Report.{ref}--{ssp}--{pop}_NUCS_msprime", zip, ref = REF, ssp = SSP, pop = POP)

#dadi_out = expand("data/dadi/RAiSD_Report.v5--{ssp}--{pop}_NUCS_msprime", zip, ssp = SSP, pop = POP)
#dadi_files = [m.replace("NUCS", n) for m in dadi_out for n in nucs]
#stats = expand("data/dadi/v5--{ssp}--{pop}_fold4_NUCS_mspms_stats.txt", zip, ssp = SSP, pop = POP)
#stats_files = [m.replace("NUCS", n) for m in stats for n in nucs]
#stats_full = expand("data/dadi/v5--{ssp}--{pop}--FULL_mspms_stats.txt", zip, ssp = SSP, pop = POP)



###########
## MUSHI ##
##########

mushi_out = expand("data/mushi/RAiSD_Report.v5--{ssp}--{pop}--msprime", zip, ssp = SSP, pop = POP)



############
## IBDseq ##
############

#data/ibdseq/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}_r2max{r2max}.hbd

#SWITCH
#ibd_files = expand("data/ibdseq/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}_r2maxR2.hbd", zip, ref = mREF, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
ibd_files = expand("data/ibdseq/v5--{ssp}--{pop}--{chrom}--{start}--{end}_r2maxR2.hbd", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
ibd_files = [f.replace("R2", r) for f in ibd_files for r in ["0.7", "0.4"]]


#################
## OTHER FILES ##
#################

        #vcf based stuff, not using
        #all_files.append(f"data/vcf/{ref}/raw/{ref}.vcf.gz")
        #all_files.appendf"data/vcf/{ref}/raw/{ref}.vcf.stats")
        #all_files.append(f"data/vcf/{ref}/filtered/{ref}_{ssp}_filtered.vcf.gz")
        #all_files.append(f"data/vcf/{ref}/filtered/{ref}_{ssp}_filtered.key")
        #all_files.append(f"data/plink/{ref}/{ref}_{ssp}_thin1M.bed")
        #all_files.append(f"data/admix/{ref}_{ssp}_thin1M.{ref_dict}.Q")
        #all_files.append(f"data/admix/{ref}_{ssp}_thin1M.{ref_dict}.P")
        #all_files.append(f"data/admix/{ref}_{ssp}_{ref_dict}_thin1M.log")


all_files = []
for ref in ["v5"]:
#for ref in ["til11", "v5"]:
    #all_files.append(f"data/trip/trip_{ref}.fa.gz")
    #all_files.append(f"data/trip/{ref}--trip--trip_nuctable.txt")
  
    all_files.append(f"data/diplo/diplo_{ref}.fa.gz")
    all_files.append(f"data/trip/trip_{ref}.fa.gz")
    all_files.append(f"data/lux/{ref}--lux--lux_nuctable.txt")
    all_files.append(f"data/diplo/{ref}--diplo--diplo_nuctable.txt")
    all_files.append(f"data/diplo/{ref}--diplo--diplo.mafs.gz")
    all_files.append(f"data/anc/{ref}_anc.fa")
    all_files.append(f"data/refs/{ref}/{ref}_FOLD")
    all_files.append(f"data/anc/{ref}_anc_FOLD")
    #all_files.append(f"data/refs/{ref}/{ref}_500M_sites.txt")
    #all_files.append(f"")
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
            #for ch in chroms:
            #    all_files.append(f"data/ngsAdmix/v5_{ssp}_{i}_{ch}_thin1M_PalmarChicoONLY.qopt")
 
#all_files.append("data/angsd_treemix/v5--trip--trip.mafscount.bed")
#print('\n'.join(pi_files))
#print(all_files)

rule all:
    input:
        ##"data/angsd_treemix/v5_treemix_filtered.fourpop.txt", #!!!
        "data/diplo/v5--diplo--diplo.mafs.gz",
        mop_files,
        all_files, 
        pi_files,
        pi_full,
        mk_files,
        relate_files,
        ibd_files,
        mushi_out,
        raisd_corrected,
        raisd_outliers,
        raisd_merged,
        pop_freqs,
        allpop_freqs
        #dadi_files,
        #stats_files,
        #stats_full,

        ##dadi_full,
        ##vcf_files,
        ##mop_final,


include: "rules/process_raw.smk"
include: "rules/mop.smk"
include: "rules/popgen.smk"
include: "rules/outgroup.smk"
include: "rules/angsd.smk"
include: "rules/mk.smk"
include: "rules/ngsRelate.smk"
include: "rules/treemix.smk"
include: "rules/angsd_vcf.smk"
include: "rules/demography.smk"
include: "rules/rdmc.smk"
