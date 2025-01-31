rule raisd:
    input:
        "data/angsd_vcf/{ref}--{ssp}--{population}--{c}--{r1}--{r2}.vcf.mop.gz"
    output:
        info = "data/raisd/RAiSD_Info.{ref}--{ssp}--{population}--{c}--{r1}--{r2}.txt",
        reports = "data/raisd/RAiSD_Report.{ref}--{ssp}--{population}--{c}--{r1}--{r2}.txt"
    params:
        population = "{ref}--{ssp}--{population}--{c}--{r1}--{r2}",
        chrom = "{c}",
        i1 = "RAiSD_Info.{ref}--{ssp}--{population}--{c}--{r1}--{r2}",
        r1 = "RAiSD_Report.{ref}--{ssp}--{population}--{c}--{r1}--{r2}.{c}"
    shell:
        """
        mkdir -p data/raisd/

        src/raisd-master/RAiSD -R -s -n {params.population} -I {input} -w 24 -m 0.05 -f 

        mv {params.i1} {output.info}
        mv {params.r1} {output.reports} 
        """

rule var_correct:
    input:
        raisd="data/raisd/RAiSD_Report.{ref}--{ssp}--{population}--{c}--{r1}--{r2}.txt",
        bed = "data/mop/{ref}--{ssp}--{population}--{c}--{r1}--{r2}.bed"
    output:
        "data/raisd/RAiSD_Report.{ref}--{ssp}--{population}--{c}--{r1}--{r2}.corrected"
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
        #awk 'BEGIN{{OFS="\\t"}}{{varnew = ($10/$9)*$5; print $1, $2, $3, $4, varnew, $6, $7, varnew*$6*$7}}' > {output}
        #awk 'BEGIN{{OFS="\\t"}}{{varnew = (1)*$5; print $1, $2, $3, $4, varnew, $6, $7, varnew*$6*$7}}' > {output}


#mushi_out = expand("data/mushi/RAiSD_Report.v5--{ssp}--{population}_msprime", zip, ssp = SSP, population = POP)
wind=50000
rule mush:
    input:
        list(set(expand("data/angsd_pi/{{ref}}--{{ssp}}--{{population}}--{chrom}--{start}--{end}.sfs", zip, chrom = mCHROM, start = mSTART, end = mEND)))
    output:
        "data/mushi/{ref}--{ssp}--{population}--mushi_ms.txt",
        "data/mushi/{ref}--{ssp}--{population}--mushi_discoal.txt"
    params:
        prefix = "data/mushi/{ref}--{ssp}--{population}--mushi",
        population_id = "{ssp}--{population}"
    #conda:
    #    "../r-environment.yml"
    shell:
        """
        python src/sim_mushi.py -i {params.population_id} -s {input} -p {params.prefix}  -n 100 -b {wind}
        """

rule demography_hbd:
    input:
        list(set(expand("data/angsd_pi/{{ref}}--{{ssp}}--{{population}}--{chrom}--{start}--{end}.sfs", zip, chrom = mCHROM, start = mSTART, end = mEND)))
    output:
        expand("data/mushi/postsims/{{ref}}--{{ssp}}--{{population}}_sim{rep}.hbd", rep = range(10))
    params:
        prefix = "data/mushi/postsims/{ref}--{ssp}--{population}",
    shell:
        "python src/sim_hbd.py -s {input} -p {params.prefix}"
        

rule mushi_raisd:
    input:
        mspms = "data/mushi/{ref}--{ssp}--{population}--mushi_ms.txt"
    output:
        info = "data/mushi/RAiSD_Info.{ref}--{ssp}--{population}--msprime",
        reports = "data/mushi/RAiSD_Report.{ref}--{ssp}--{population}--msprime"
    params:
        population = "{ref}--{ssp}--{population}--msprime",
        i1 = "RAiSD_Info.{ref}--{ssp}--{population}--msprime",
        r1 = "RAiSD_Report.{ref}--{ssp}--{population}--msprime"
    shell:
        """
        src/raisd-master/RAiSD -I {input.mspms} -n {params.population} -L {wind} -w 24 -m 0.05 
        mv {params.i1} {output.info}
        mv {params.r1} {output.reports} 
        """

rule raisd_outliers:
    input:
        msprime = "data/mushi/RAiSD_Report.{ref}--{ssp}--{population}--msprime",
        raisd = "data/raisd/RAiSD_Report.{ref}--{ssp}--{population}--{c}--{r1}--{r2}.corrected"
    output:
        "data/raisd/RAiSD_Report.{ref}--{ssp}--{population}--{c}--{r1}--{r2}.corrected_block_outliers"
    shell:
        """
        quantile=`grep -v "//" {input.msprime} | awk '{{print $2}}' | sort -rg | perl -e '$d=0.05;@l=<>;print $l[int($d*$#l)]'`
        cat {input.raisd} | awk -v quantile=$quantile '{{OFS = "\\t"}}; $8 > quantile {{print $1, $2-1, $2, $5, $6, $7, $8}}' | bedtools sort -i stdin  | bedtools merge -i stdin -d 100000 -c 7 -o max > {output}
        """

rule sweep_regions:
    input:
        mushi = "data/mushi/RAiSD_Report.{ref}--{ssp}--{population}--msprime",
        raisd = "data/raisd/RAiSD_Report.{ref}--{ssp}--{population}--{c}--{r1}--{r2}.corrected"
    output:
        beds = "data/sweep_regions/sweeps.{ref}--{ssp}--{population}--{c}--{r1}--{r2}.bed",
         outs = "data/sweep_regions/sweeps.{ref}--{ssp}--{population}--{c}--{r1}--{r2}.csv"
    params:
         cutoff = 1 - 0.001,
         merge = 50000
    shell:
        """
        #get outlier from sims
        cutoff=`src/sweep_regions/sweep_regions.R neutral_cutoff \
          --quantile {params.cutoff} \
          -f {input.mushi} \
          --input_raisd2 TRUE \
          --delimeter "\\t" | cut -d ":" -f2`

        #run sweep_regions
        src/sweep_regions/sweep_regions.R sweep_regions \
        --data_frame {input.raisd} \
        -R TRUE \
        --delimeter "\\t" \
        --merge_size {params.merge} \
        --cutoff ${{cutoff}} \
        --min_size 1e3 \
        --max_size 1e8 \
        --out_file {output.outs}
        """


no_pc_bed = [i.replace(".csv", ".bed") for i in population_sweep_regions]
rule shared:
    input:
        no_pc_bed
    output:
        "data/sweep_regions/{ref}--allpopulations--shared_outliers.txt"
    resources:
        mem_mb = 1000
    shell:
        """
        cat {input} |\
        bedtools sort -i stdin |\
        bedtools intersect -a stdin -b {input} -filenames -wb |\
        cut -f1-4 |\
        bedtools sort -i stdin |\
        bedtools merge -c 4 -o distinct -delim "," | awk '{{print $1 "\\t" $2 "\\t" $3 "\\t" $3-$2 "\\t" $4}}' > {output} 
        """
        
#rule merge_outliers:
#    input:
#        allfiles = list(set(expand("data/raisd/RAiSD_Report.{{ref}}--{{ssp}}--{{population}}--{c}--{r1}--{r2}.corrected_block_outliers", zip, c = mCHROM, r1 = mSTART, r2 = mEND)))
#    output:
#        "data/raisd/{ref}--{ssp}--{population}.corrected_block_outliers_merged.txt"
#    shell:
#        "cat {input.allfiles} > {output}"
     
#files that should be considered for merging, no full palmar chicos 
#populations_string = "data/raisd/v5--LR--Amatlan_de_Canas.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--LR--Crucero_Lagunitas.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--LR--Los_Guajes.corrected_block_outliers_merged.txt  " + \
#"data/raisd/v5--LR--random1_Palmar_Chico.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--LR--random2_Palmar_Chico.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--LR--random.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--LR--San_Lorenzo.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--Teo--Amatlan_de_Canas.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--Teo--Crucero_Lagunitas.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--Teo--El_Rodeo.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--Teo--Los_Guajes.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--Teo--random1_Palmar_Chico.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--Teo--random2_Palmar_Chico.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--Teo--random.corrected_block_outliers_merged.txt " + \
#"data/raisd/v5--Teo--San_Lorenzo.corrected_block_outliers_merged.txt"

#expand("data/raisd/{{ref}}--{ssp}--{population}.corrected_block_outliers_merged.txt", zip, ssp  = mSSP, population = mPOP)
#rule shared:
#    input:
#        raisd_merged
#    output:
#        "data/raisd/{ref}--allpopulations--shared_outliers.txt"
#    shell:        
#        """
#        cat {input} |\
#        bedtools sort -i stdin |\
#        bedtools intersect -a stdin -b {input} -filenames -wb |\
#        cut -f1-5 |\
#        bedtools sort -i stdin |\
#        bedtools merge -c 5 -o distinct -delim "," | awk '{{print $1 "\\t" $2 "\\t" $3 "\\t" $3-$2 "\\t" $4}}' > {output} 
#        """


#cat {populations_string} |\
#bedtools sort -i stdin |\
#bedtools merge -i stdin |\
#bedtools intersect -a stdin -b {populations_string} -filenames -wb |\
#cut -f1-4 |\
#bedtools sort -i stdin |\
#bedtools merge -c 4 -o distinct -delim "," | awk '{{print $1 "\\t" $2 "\\t" $3 "\\t" $3-$2 "\\t" $4}}' > {output} 



