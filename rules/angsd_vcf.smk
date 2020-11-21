#vcf per pop
rule raisd_beagle:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        bams = "data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt"
    output:
        "data/angsd_vcf/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.vcf.gz"
    params:
        scratch = my_scratch,
        prefix = my_scratch + "{ref}--{ssp}--{pop}--{c}--{r1}--{r2}",
        final = "data/angsd_vcf/",
        chrom = "{c}",
        r1 = "{r1}",
        r2 = "{r2}"
    shell:
        """
        mkdir -p {params.scratch}
        module load bio
        angsd -b {input.bams} -P 5 \
        -ref {input.ref} \
        -r {params.chrom}:{params.r1}-{params.r2} \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 30 -minQ 30 \
        -doCounts 1 -setMinDepth 5 -setMaxDepth 150 -doGeno 3 -dovcf 1 -gl 2 -dopost 2 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -out {params.prefix}
        
        mv {params.prefix}.vcf.gz {params.final}
        mv {params.prefix}.mafs.gz {params.final}
        mv {params.prefix}.arg {params.final}
        """

#data/mop/v5--LR--random1_Palmar_Chico--chr1--0--308452471.bed
rule vcf_mop:
    input:
        mop = "data/mop/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.bed",
        vcf = "data/angsd_vcf/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.vcf.gz"
    output:
        "data/angsd_vcf/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}_mop.vcf"
    shell:
        "bedtools intersect -header -a {input.vcf} -b {input.mop} > {output}"

rule IBDseq:
    input:
        "data/angsd_vcf/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}_mop.vcf"
    output:
        "data/ibdseq/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}_r2max{r2max}.hbd"
    params:
        prefix = "data/ibdseq/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}_r2max{r2max}",
        r2max = "{r2max}"
    shell:
        "java -jar src/ibdseq.r1206.jar gt={input} out={params.prefix} nthreads=12 r2max={params.r2max}"

