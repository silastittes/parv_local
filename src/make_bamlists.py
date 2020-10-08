from glob import glob
import pandas as pd

pop_key = pd.read_csv("pop_key", sep = "\t")

chrom_key = pd.read_csv("chroms.txt", sep = "\t")

pops = list(set(pop_key['geo']))
ssp = list(set(pop_key['species']))

paths = {"til11" : "data/bams/til11-alignments/deduped-bam/", "v5": "data/bams/v5-alignments/deduped-bam/"}

for pop in pops:
    for sp in ssp:
        if pop == "El_Rodeo" and sp == "LR":
            pass
        else:
            for ref, path in paths.items():
                plist = list(pop_key.query(f"geo == '{pop}' & species == '{sp}'")["JRIAL_ID"])
                filenames = [f"{path}{ind}.deduped.bam" for ind in plist]
                outfile = open(f"data/bamlist/{ref}--{sp}--{pop}__bamlist.txt", "w")
                popfile = open(f"data/bamlist/{ref}--{sp}--{pop}__ID.txt", "w")
                print('\n'.join(filenames), file = outfile)
                print('\n'.join(plist), file = popfile)

            for index, row in chrom_key.iterrows():
                file_str = open(f"data/mop/{row['ref']}--{sp}--{pop}--{row['chrom']}--{row['start']}--{row['end']}.txt", "w")
                print("Oh hello. The file name is all that matters here.", file = file_str)

#make names that indicate the ref, pop, chromosome  combinations for mop!!!!


#lr_bam = glob("data/B73v5-alignments/deduped-bam/*bam")
#teo_bam = glob("data/B73v5-alignments/deduped-bam/*bam")
