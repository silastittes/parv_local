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
        sites = til_raw_stats,
        genos = til_raw_stats + "genotypes"
    shell:
        "tidy_vcf -t 50000 -v {input} -o {output.sites} -g {output.genos}"


rule filtered_stat:
    input:
        til_filtered
    output:
        sites = til_filtered_stats,
        genos = til_filtered_stats + "genotypes"
    shell:
        "tidy_vcf -t 50000 -v {input} -o {output.sites} -g {output.genos}"



rule filter:
    input:
        til_raw
    output:
        til_filtered
    shell:
        """
        bcftools filter -O z -o {output} -e 'COUNT(GT=="mis") >= 0.2*195 | QUAL < 30 | MIN(FORMAT/DP[@{teo_key}]) < 10 |  MAX(N_ALT) > 1 | TYPE!="snp" | ExcessHet > 20 | -0.5 > InbreedingCoeff > 0.5' {input}
        """
