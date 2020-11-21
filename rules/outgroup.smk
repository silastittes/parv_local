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


