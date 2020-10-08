rule admixture:
    input:
        key = "data/vcf/{ref}/filtered/{ref}_{ssp}_filtered.key",
        bed = "data/plink/{ref}/{ref}_{ssp}_thin1M.bed"
    output:
        P = "data/admix/{ref}_{ssp}_thin1M.{K}.P",
        Q = "data/admix/{ref}_{ssp}_thin1M.{K}.Q",
        logfile = "data/admix/{ref}_{ssp}_{K}_thin1M.log",
    params:
        K = "{K}",
        ssp = "{ssp}",
        ref = "{ref}"
    shell:
        """
        admixture {input.bed} {params.K} > {output.logfile}
        mv {params.ref}_{params.ssp}_thin1M.{params.K}.Q {output.Q}
        mv {params.ref}_{params.ssp}_thin1M.{params.K}.P {output.P}
        """ 

#data/bamlist/til11--Teo--Los_Guajes__ID.txt
rule raw_pi:
    input:
        vcf = "data/vcf/{ref}/filtered/{ref}_{ssp}_filtered.vcf.gz",
    output:
        "data/pi/{ref}--{ssp}--{pop}__pi.txt"
    params:
        ref = "{ref}",
        ssp = "{ssp}",
        pop = "{pop}"
    shell:
        "vcftools --keep data/bamlist/{params.ref}--{params.ssp}--{params.pop}__ID.txt --gzvcf {input.vcf} --window-pi 10000 --stdout > {output}"


