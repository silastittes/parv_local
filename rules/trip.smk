rule trip_bam:
    input:
        r1 = "data/trip/Tripsacum_KD17021399_HHYHYBCXY_L2_1.fq.gz",
        r2 = "data/trip/Tripsacum_KD17021399_HHYHYBCXY_L2_2.fq.gz",
        ref = "data/refs/{ref}/{ref}.fa"
    output:
        bam = temp("data/trip/trip_{ref}.bam")
    shell:
        """
        bwa mem -t 48 {input.ref} {input.r1} {input.r2} > {output.bam}
        """

rule trip_sort:
    input:
        "data/trip/trip_{ref}.bam",
    output:
        "data/trip/trip_{ref}.sorted.bam",
    shell:
        """
        samtools sort --threads 48 -o {output} {input}
        """

rule trip_fasta:
    input:
        "data/trip/trip_{ref}.sorted.bam",
    output:
        "data/trip/trip_{ref}.fa.gz",
    params:
        prefix = "data/trip/trip_{ref}",
    shell:
        """
        module load angsd
        angsd -i {input} -doFasta 2 -doCounts 1 -out {params.prefix}
        """

rule trip_beagle:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        trip = "data/trip/trip_{ref}.fa.gz",
        bams = "data/trip/{ref}--{ssp}--{pop}__bamlist.txt"
    output:
        mafs = "data/angsd_pi/{ref}--{ssp}--{pop}.mafs.gz"
    params:
        prefix = "data/angsd_pi/{ref}--{ssp}--{pop}"
    shell:
        """
        module load angsd
        angsd -GL 1 -P 5 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 30 -minQ 30 \
        -ref {input.ref}  -anc {input.trip} \
        -doMaf 2 -doMajorMinor 4 \
        -bam {input.bams} -out {params.prefix}
        """

