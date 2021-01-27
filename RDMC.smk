#separate snakemake file to run rdmc
import pandas as pd
import numpy as np


CHROM, START, END = glob_wildcards("data/rdmc/freq/v5--allpops--{chrom}--{start}--{end}_freq.txt.gz")

chrom_dict = {}
for c, s, e in zip(CHROM, START, END):
    chrom_dict[c] = f"{c}--{s}--{e}"

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



sweep_df = pd.read_csv("data/raisd/v5--allpops--shared_outliers.txt",
            sep = "\t",
            names = ["chrom", "start", "end", "length", "files"])

sweep_idx = []
for sweep in sweep_df['files']:
    sweep_files = sweep.split(',')
    s_idx = [f"{i+1}" for i,x in enumerate(FREQ_POPS) for s in sweep_files if s.find(x) != -1]
    if len(s_idx) < 2 or len(s_idx) > len(FREQ_POPS)-2:
        sweep_idx.append("NA")
    else:
        sweep_idx.append('-'.join(s_idx))

sweep_df['sweep_idx'] = sweep_idx

out_file = [f"data/rdmc/fitted/v5--sweep_{chrom_dict[row['chrom']]}_start{str(row['start'])}_end{str(row['end'])}_pops{row['sweep_idx']}.txt" for index, row in sweep_df.iterrows()]
sweep_df['out_file'] = out_file
sweep_df = sweep_df[sweep_df['sweep_idx'] != "NA"].reset_index() 

neutrals = [f"data/rdmc/freq/v5--NEUTRAL--{chrom}--{start}--{end}_freq.txt.gz" for chrom, start, end in list(set(zip(CHROM, START, END)))]

print('\n'.join(neutrals))


###########
## RULES ##
###########

rule all:
    input:
        list(sweep_df['out_file']),
        "data/rdmc/v5--neutral_freqs.txt"
        #neutrals

rule get_sweep:
    input:
        "data/rdmc/freq/{ref}--allpops--{chrom}--{start}--{end}_freq.txt.gz"
    output:
        "data/rdmc/sweep_freq/{ref}--sweep--{chrom}--{start}--{end}_start{sweep_start}_end{sweep_end}_pops{pops}.txt"
    params:
        chrom = "{chrom}",
        start = "{sweep_start}",
        end = "{sweep_end}"
    shell:
        """
        zcat {input} |  awk -v chrom={params.chrom} -v start={params.start} -v end={params.end} '$1==chrom && $2 >= start && $3 <= end {{print $0}}' > {output}
        """


rule get_neutral:
    input:
        expand("data/rdmc/freq/{{ref}}--allpops--{chrom}--{start}--{end}_freq.txt.gz", zip, chrom = CHROM, start=START, end=END)
    output:
        "data/rdmc/{ref}--neutral_freqs.txt"
    shell:
        "> {output}; for i in `echo {input}`; do zcat $i | shuf -n 10000 >> {output}; done"
    

#"sweep_chr1--0--308452471_start995115_end1294915_pops3-4-5-6-8-9-10-11-12.txt"
rule rdmc_cli:
    input:
        gmap = "data/map/ogut_{ref}.map.txt",
        neutral_file = "data/rdmc/{ref}--neutral_freqs.txt",
        sweep_file = "data/rdmc/sweep_freq/{ref}--sweep--{chrom}--{start}--{end}_start{sweep_start}_end{sweep_end}_pops{pops}.txt",
    output:
        "data/rdmc/fitted/{ref}--sweep_{chrom}--{start}--{end}_start{sweep_start}_end{sweep_end}_pops{pops}.txt"
    params:
        pops = "{pops}",
        start = "{sweep_start}",
        end = "{sweep_end}"
    conda:
        "rdmc-environment.yml"
    shell:
        """
        Rscript src/rdmc_cli.R --neutral_file {input.neutral_file} --sweep_file {input.sweep_file} --pop_ids {params.pops} --gmap {input.gmap} --out_file {output}
        """

