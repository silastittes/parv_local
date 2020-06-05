import re

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


REF, REGIONS = glob_wildcards("data/vcf/{ref}/vcffiles_split_raw/{regions}.vcf.gz")

REGIONS = natural_sort(REGIONS)

til_raw = "data/vcf/til11/raw/til11.vcf.gz"
til_raw_stats = til_raw.replace("gz", "stats")

rule all:
    input:
        til_raw, til_raw_stats


include: "rules/process_raw.smk"
