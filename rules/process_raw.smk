
#til_raw = "data/vcf/til11/raw/til11.vcf.gz"
#til_raw_stats = til_raw.replace("gz", "stats")
#til_filtered = "data/vcf/til11/filtered/til11_filtered.vcf.gz"

#rule vcf_concat:
#    input:
#        #expand("data/vcf/{ref}/vcffiles_split_raw/{regions}.vcf.gz", zip, regions = REGIONS, ref = REF)
#        i=expand("data/vcf/{ref}/vcffiles_split_raw/{regions}.vcf", zip, regions = REGIONS, ref = REF)
#    output:
#        "data/vcf/{ref}/raw/{ref}.vcf.gz"
#    shell:
#        "bcftools concat {input.i} -O z -o {output}"
#        #"bcftools concat data/vcf/til11/vcffiles_split_raw/*.vcf -O z -o {output}"

rule raw_stat:
    input:
        "data/vcf/{ref}/raw/{ref}.vcf.gz"
    output:
        sites = "data/vcf/{ref}/raw/{ref}.vcf.stats",
        genos = "data/vcf/{ref}/raw/{ref}.vcf.genotypes"
    shell:
        "tidy_vcf -t 10000 -v {input} -o {output.sites} -g {output.genos}"


#need to update to allow different pop keys
rule filter:
    input:
        vcf = "data/vcf/{ref}/raw/{ref}.vcf.gz",
    output:
        "data/vcf/{ref}/filtered/{ref}_{ssp}_filtered.vcf.gz"
    params:
        key = "{ssp}"
        
    shell:
        """
        vcftools --keep {params.key}_key.txt --gzvcf {input.vcf} --max-alleles 2 --min-alleles 2 --remove-indels --minDP 5 --maxDP 100 --minGQ 30 --minQ 30 --max-missing 0.7 --recode --stdout | gzip > {output}
        """

rule IDs:
    input:
        "data/vcf/{ref}/filtered/{ref}_{ssp}_filtered.vcf.gz"
    output:
        "data/vcf/{ref}/filtered/{ref}_{ssp}_filtered.key"
    shell:
        "bcftools query -l {input} > {output}"
        

rule thin_plink:
    input:
        "data/vcf/{ref}/filtered/{ref}_{ssp}_filtered.vcf.gz"
    params:
        prefix = "data/plink/{ref}/{ref}_{ssp}_thin1M"
    output:
        "data/plink/{ref}/{ref}_{ssp}_thin1M.bed"
    shell:
        """
        module load plink/1.90
        plink --make-bed --thin-count 1000000 --vcf {input}  --out {params.prefix}
        """

#data/trip/trip_v5.fa.gz  data/trip/trip_v5_FOLD
rule fold:
    input:
        ref = "data/trip/trip_{ref}.fa",
        #ref = "data/refs/{ref}/{ref}.fa",
        gff = "data/refs/{ref}/{ref}.gff3.gz"
    output:
        fold = "data/trip/trip_{ref}_FOLD",
        fold_log = "data/trip/trip_{ref}_FOLD.log"
    shell:
        """
        python src/cds_fold/cds_fold.py {input.ref} {input.gff} -o {output.fold} > {output.fold_log}
        awk '$3 ~ /4/' {output.fold} > {output.fold}4
        awk '$3 ~ /0/' {output.fold} > {output.fold}0
        """ 
            
