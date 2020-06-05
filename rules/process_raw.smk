rule vcf_concat:
    input:
        expand("data/vcf/{ref}/vcffiles_split_raw/{regions}.vcf.gz", zip, regions = REGIONS, ref = REF)
    output:
        til_raw
    shell:
        "bcftools concat {input} -O z -o {output}"

rule raw_stat:
    input:
        til_raw
    output:
        til_raw_stats
    shell:
        "bcftools stats -s - {input} > {output}"
