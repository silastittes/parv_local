#vcf per pop
rule raisd_beagle:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        bams = "data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt"
    output:
        "data/angsd_vcf/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.vcf.gz"
    params:
        prefix = "data/angsd_vcf/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}",
        chrom = "{c}",
        r1 = "{r1}",
        r2 = "{r2}"
    shell:
        """
        module load bio
        angsd -b {input.bams} -P 5 \
        -ref {input.ref} \
        -r {params.chrom}:{params.r1}-{params.r2} \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 30 -minQ 30 \
        -dovcf 1 -gl 1 -dopost 1 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -out {params.prefix}
        """

rule raisd:
    input:
        "data/angsd_vcf/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.vcf.gz"
    output:
        i1 = "RAiSD_Info.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}",
        r1 = "RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}",
        info = "data/raisd/RAiSD_Info.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}",
        report = "data/raisd/RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}"
    params:
        pop = "{ref}--{ssp}--{pop}--{c}--{r1}--{r2}"
    shell:
        """
        mkdir -p data/raisd/

        src/raisd-master/RAiSD \
          -R -s \
          -m 0.05 \
          -n {params.pop} \
          -I {input}

        mv {output.i1} {output.info}
        mv {output.r1} {output.report} 
        """


rule var_correct:
    input:
        raisd="data/raisd/RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}",
        bed = "data/mop/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.bed"
    output:
        "data/raisd/RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.corrected"
    params:
        c = "{c}"
    shell:
        """
        awk -v c={params.c} 'BEGIN{{OFS = "\\t"}}; NR > 1 {{print c, $2-1, $3, $0, $3-$2+1}}' {input.raisd} |\
        bedtools groupby -i - -g 1,2,3,4,5,6 -c 7,8,9,10,11 -o last,last,last,last,distinct |\
        bedtools intersect -a - -b {input.bed} |\
        awk 'BEGIN{{OFS = "\\t"}};{{print $0, $3-$2}}' |\
        bedtools groupby -i stdin -g 1,4,5,6 -c 7,8,9,10,11,12 -o distinct,distinct,distinct,distinct,distinct,sum |\
        awk 'BEGIN{{OFS="\\t"}}{{varnew = ($10/$9)*$5; print $1, $2, $3, $4, varnew, $6, $7, varnew*$6*$7}}' > {output}
        """
