#data/trip/trip_v5.fa.gz

rule trip_bam:
    input:
        r1 = "data/{out}/{out}_1.fastq.gz",
        r2 = "data/{out}/{out}_2.fastq.gz",
        ref = "data/refs/{ref}/{ref}.fa"
    output:
        bam = temp("data/{out}/{out}_{ref}.bam")
    shell:
        """
        bwa mem -t 48 {input.ref} {input.r1} {input.r2} > {output.bam}
        """

rule trip_sort:
    input:
        "data/{out}/{out}_{ref}.bam",
    output:
        "data/{out}/{out}_{ref}.sorted.bam",
    shell:
        """
        samtools sort --threads 48 -o {output} {input}
        """

rule trip_fasta:
    input:
        "data/{out}/{out}_{ref}.sorted.bam",
    output:
        "data/{out}/{out}_{ref}.fa.gz",
    params:
        prefix = "data/{out}/{out}_{ref}",
    shell:
        """
        module load angsd
        angsd -i {input} -doFasta 2 -doCounts 1 -out {params.prefix}
        """

rule trip_beagle:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        trip = "data/{out}/{out}_{ref}.fa.gz",
        bams = "data/{out}/{ref}--{ssp}--{pop}__bamlist.txt"
    output:
        mafs = "data/{out}/{ref}--{ssp}--{pop}.mafs.gz"
    params:
        prefix = "data/{out}/{ref}--{ssp}--{pop}"
    shell:
        """
        module load angsd
        angsd -GL 1 -P 5 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 30 -minQ 30 \
        -ref {input.ref} \
        -doCounts 1 -doMaf 2 -doMajorMinor 4 \
        -bam {input.bams} -out {params.prefix}
        """
        #-setMinDepthInd 5 \

rule nuctable:
    input:
       "data/{out}/{ref}--{ssp}--{pop}.mafs.gz" 
    output:
        "data/{out}/{ref}--{ssp}--{pop}_nuctable.txt"
    shell:
        "python src/maf2nuccounts.py -f {input}  -t unknown -m 0.001  > {output}"


grps = ["lux", "diplo"]
rule anc_bed:
    input:
        #"data/diplo/{ref}--diplo--diplo_nuctable.txt"
        expand("data/{out}/{{ref}}--{ssp}--{pop}_nuctable.txt", zip, out = grps, ssp = grps, pop = grps)
    output:
        "data/anc/{ref}_anc.bed"
    shell:
        """
        bedtools unionbedg -filler MISSING -i {input} | awk '$4 !~ /,/ && $5 !~ /,/ {{print $0}}' | awk '{{if($4 == $5){{print $1 "\\t" $2 "\\t" $3 "\\t" $4}}; if($4 ~ /MISSING/){{print $1 "\\t" $2 "\\t" $3 "\\t" $5}}}}' > {output}
        """

#bedtools unionbedg -filler MISSING -i {input} | awk '$4 !~ /,/ && $5 !~ /,/ {{print $0}}' | awk '{{if($4 == $5){{print $1 "\\t" $2 "\\t" $3 "\\t" $4}}; if($4 ~ /MISSING/){{print $1 "\\t" $2 "\\t" $3 "\\t" $5}}; if($5 ~/MISSING/){{print $1 "\\t" $2 "\\t" $3 "\\t" $4}}}}' > {output}
#alternatives for getting ancestral nucleotides
#awk '$4 !~ /,/{{print $0}}' {input} > {output}    
#bedtools unionbedg -filler 'MISSING' -i {input} | grep -v MISSING | awk '$4 ~ /^[ATGC]$/ && $5 ~ /^[ATGC]$/ && $4 == $5 {{print $1 "\\t" $2 "\\t" $3 "\\t" $4}}' > {output}

rule anc_ref:
    input:
        gbed = "data/refs/{ref}/{ref}.gbed",
        bed = "data/anc/{ref}_anc.bed"
    output:
        "data/anc/{ref}_anc.fa"
    shell:
        "python src/impute_fasta.py -g {input.gbed} -b {input.bed} > {output}"


rule fold_anc:
    input:
        ref = "data/anc/{ref}_anc.fa",
        gff = "data/refs/{ref}/{ref}.gff3.gz"
    output:
        fold = "data/anc/{ref}_anc_FOLD",
        fold4 = "data/anc/{ref}_anc_FOLD4",
        fold0 = "data/anc/{ref}_anc_FOLD0",
        fold_log = "data/anc/{ref}_anc_FOLD.log"
    shell:
        """
        python src/cds_fold/cds_fold.py -c {input.ref} {input.gff} -o {output.fold} > {output.fold_log}
        awk '$3 ~ /4/' {output.fold} > {output.fold}4
        awk '$3 ~ /0/' {output.fold} > {output.fold}0
        """




rule diplomaf2count:
    input:
        "data/diplo/{ref}--diplo--diplo.mafs.gz"
    output:
        "data/angsd_treemix/{ref}--diplo--diplo.mafscount.bed"
    shell:
        """
        zcat < {input} | awk 'NR>1{{scale=1/(2*$8); freq = $7+scale/2 -($7+scale/2)%scale; print $1 "\\t" $2-1 "\\t" $2 "\\t" int(freq*$8*2) "," int((1-freq)*$8*2)}}' > {output}
        """



#        mkdir -p {params.scratch}
#        module load bio
#        angsd -b {input.bams} -P 5 \
#        -ref {input.ref} \
#        -r {params.chrom}:{params.r1}-{params.r2} \
#        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 30 -minQ 30 \
#        -doCounts 1 -setMinDepth 5 -setMaxDepth 150 -doGeno 3 -dovcf 1 -gl 2 -dopost 2 -domajorminor 1 -domaf 1 -snp_pval 1e-6 -out {params.prefix}
        
 #       mv {params.prefix}.vcf.gz {params.final}
 #       mv {params.prefix}.mafs.gz {params.final}
 #       mv {params.prefix}.arg {params.final}


