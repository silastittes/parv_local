import re

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


REF, REGIONS = glob_wildcards("data/vcf/{ref}/vcffiles_split_raw/{regions}.vcf")

REGIONS = natural_sort(REGIONS)

K_dict = {"Teo_k": 6, "LR_k": 5}

#######
# mop #
#######
#data/bamlist/til11--LR--Amatlan_de_Canas__bamlist.txt
REF, SSP, POP = glob_wildcards("data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt")

mREF, mSSP, mPOP, mCHROM, mSTART, mEND = glob_wildcards("data/mop/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.txt")

mop_files = expand("data/mop/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.bed", zip, 
                    ref = mREF, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)

mop_final = expand("data/mop/{ref}--{ssp}--{pop}.bed", zip, ref = mREF, ssp = mSSP, pop = mPOP)


##########
##  PI  ##
##########

#pi_files = expand("data/angsd_pi/{ref}--{ssp}--{pop}.{win}kb_theta.thetasWindow.gz.pestPG", zip, ref = REF, pop = POP, ssp = SSP)

pi_files = expand("data/angsd_pi/v5--{ssp}--{pop}.WINDOWBP_theta.thetasWindow.gz.pestPG", zip, ssp = SSP, pop = POP)

pi_files = [p.replace("WINDOW", w) for p in pi_files for w in ["1000", "100000", "1000000"]]



##########
### MK ###
##########


#0: A,T -> G,C
#0: G,C -> A,T
#0: A,T -> A,T + G,C -> G,C

fold = ["0", "4"]
nucs = ["AG,AC,TG,TC",  "GA,GT,CA,CT", "AA,AT,TA,TT,GG,GC,CG,CC"]
#mk_sites = expand("data/angsd_sfs/{ref}--{ssp}--{pop}_foldFOLD_NUCS_sfs.txt",  zip, ref = REF, ssp = SSP, pop = POP)
mk_sites = expand("data/angsd_sfs/v5--{ssp}--{pop}_foldFOLD_NUCS_sfs.txt",  zip, ssp = SSP, pop = POP)

mk_files = [m.replace("FOLD", f).replace("NUCS", n) for m in mk_sites for f in fold for n in nucs]

print('\n'.join(mk_files))

all_files = []
#for ref in ["til11"]:
for ref in ["v5"]:
#for ref in ["til11", "v5"]:
    for ssp in ["Teo", "LR"]:
        #ref_dict = K_dict[f'{ssp}_k']
        #all_files.append(f"data/vcf/{ref}/raw/{ref}.vcf.gz")
        #all_files.appendf"data/vcf/{ref}/raw/{ref}.vcf.stats")
        #all_files.append(f"data/vcf/{ref}/filtered/{ref}_{ssp}_filtered.vcf.gz")
        #all_files.append(f"data/vcf/{ref}/filtered/{ref}_{ssp}_filtered.key")
        #all_files.append(f"data/plink/{ref}/{ref}_{ssp}_thin1M.bed")
        #all_files.append(f"data/admix/{ref}_{ssp}_thin1M.{ref_dict}.Q")
        #all_files.append(f"data/admix/{ref}_{ssp}_thin1M.{ref_dict}.P")
        #all_files.append(f"data/admix/{ref}_{ssp}_{ref_dict}_thin1M.log")

        #!!!
        all_files.append(f"data/beagle/{ref}--{ssp}.beagle.gz")
    all_files.append(f"data/refs/{ref}/{ref}_FOLD")
    all_files.append(f"data/trip/trip_{ref}_FOLD")


#print('\n'.join(pi_files))
print(all_files)

rule all:
    input:
        #"data/trip/trip_til11.fa.gz",
        "data/trip/trip_v5.fa.gz",
        "data/refs/v5/v5_FOLD",
        "data/trip/trip_v5_FOLD",
        pi_files,
        mk_files #,all_files
        #mop_files,
        #mop_final,


include: "rules/process_raw.smk"
include: "rules/mop.smk"
include: "rules/popgen.smk"
include: "rules/trip.smk"
include: "rules/angsd.smk"
include: "rules/mk.smk"
