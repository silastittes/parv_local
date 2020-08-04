import re

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

teo_key = "teo.txt" 
lr_key = "lr.txt"

REF, REGIONS = glob_wildcards("data/vcf/{ref}/vcffiles_split_raw/{regions}.vcf.gz")

REGIONS = natural_sort(REGIONS)

til_raw = "data/vcf/til11/raw/til11.vcf.gz"
til_raw_stats = til_raw.replace("gz", "stats")
til_filtered = "data/vcf/til11/filtered/til11_filtered.vcf.gz"
til_filtered_stats = "data/vcf/til11/filtered/til11_filtered.vcf.gz".replace("gz", "stats")


rule all:
    input:
        til_raw, til_raw_stats,
        til_filtered, til_filtered_stats


include: "rules/process_raw.smk"
