rule raisd:
    input:
        "data/angsd_vcf/{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.vcf.gz"
    output:
        info = "data/raisd/RAiSD_Info.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.txt",
        report = "data/raisd/RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.txt"
    params:
        pop = "{ref}--{ssp}--{pop}--{c}--{r1}--{r2}",
        chrom = "{c}",
        i1 = "RAiSD_Info.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}",
        r1 = "RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.{c}"
    shell:
        """
        mkdir -p data/raisd/

        src/raisd-master/RAiSD \
          -R -s \
          -m 0.05 \
          -n {params.pop} \
          -I {input}

        mv {params.i1} {output.info}
        mv {params.r1} {output.report} 
        """


rule var_correct:
    input:
        raisd="data/raisd/RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.txt",
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


rule dadi:
    input:
        "data/angsd_sfs/{ref}--{ssp}--{pop}_fold4_{nucs}_sfs.txt"
    output: 
        "data/dadi/{ref}--{ssp}--{pop}_fold4_{nucs}_msprime.txt"
    params:
        prefix = "data/dadi/{ref}--{ssp}--{pop}_fold4_{nucs}"
    shell:
        "python src/sim_dadi_raisd.py -n 1000 -b 10000 -s {input} -p {params.prefix}"

rule sample_stats:
    input:
        "data/dadi/{ref}--{ssp}--{pop}_fold4_{nucs}_msprime.txt"
    output:
        "data/dadi/{ref}--{ssp}--{pop}_fold4_{nucs}_mspms_stats.txt"
    shell:
        "cat data/dadi/v5--LR--random1_Palmar_Chico_fold4_GA,GT,CA,CT_msprime.txt| src/msdir/sample_stats > {output}"

rule msprime_raisd:
    input:
        mspms = "data/dadi/{ref}--{ssp}--{pop}_fold4_{nucs}_msprime.txt",
        stats = "data/dadi/{ref}--{ssp}--{pop}_fold4_{nucs}_mspms_stats.txt"
    output:
        info = "data/dadi/RAiSD_Info.{ref}--{ssp}--{pop}_{nucs}_msprime",
        report = "data/dadi/RAiSD_Report.{ref}--{ssp}--{pop}_{nucs}_msprime"
    params:
        pop = "{ref}--{ssp}--{pop}_{nucs}_msprime",
        i1 = "RAiSD_Info.{ref}--{ssp}--{pop}_{nucs}_msprime",
        r1 = "RAiSD_Report.{ref}--{ssp}--{pop}_{nucs}_msprime"
    shell:
        """
        src/raisd-master/RAiSD -I {input.mspms} -n {params.pop} -L 10000
        mv {params.i1} {output.info}
        mv {params.r1} {output.report} 
        """

rule raisd_outliers:
    input:
        msprime = "data/dadi/RAiSD_Report.{ref}--{ssp}--{pop}_AT,TA,GC,CG_msprime",
        raisd = "data/raisd/RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.corrected"
    output:
        "data/raisd/RAiSD_Report.{ref}--{ssp}--{pop}--{c}--{r1}--{r2}.corrected_block_outliers"
    shell:
        """
        quantile=`grep -v "//" notebooks/RAiSD_Report.test | awk '{{print $2}}' | sort -rg | perl -e '$d=0.01;@l=<>;print $l[int($d*$#l)]'`
        cat {input} | awk -v quantile=$quantile '{{OFS = "\\t"}}; $8 > quantile {{print $1, $2-1, $2, $5, $6, $7, $8}}' | bedtools sort -i stdin  | bedtools merge -i stdin -d 100000 -c 7 -o max > {output}
        """

rule merge_outliers:
    input:
        allfiles = list(set(expand("data/raisd/RAiSD_Report.{{ref}}--{{ssp}}--{{pop}}--{c}--{r1}--{r2}.corrected_block_outliers", zip, c = mCHROM, r1 = mSTART, r2 = mEND)))
    output:
        "data/raisd/{ref}--{ssp}--{pop}.corrected_block_outliers_merged.txt"
    shell:
        "cat {input.allfiles} > {output}"
