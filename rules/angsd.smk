
rule full_beagle:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        trip = "data/trip/trip_{ref}.fa.gz",
        bams = "data/bamlist/{ref}-{ssp}-allBams.txt"
    output:
        "data/beagle/{ref}--{ssp}.beagle.gz"
    params:
        prefix = "data/beagle/{ref}--{ssp}"
    shell:
        """
        module load angsd
        angsd -GL 1 -P 10 \
        -ref {input.ref}  -anc {input.trip} \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 30 -minQ 30 \
        -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam {input.bams} -out {params.prefix}
        """

rule pop_beagle:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        trip = "data/trip/trip_{ref}.fa.gz",
        bams = "data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt"
    output:
        saf = "data/angsd_pi/{ref}--{ssp}--{pop}.saf.gz",
    params:
        prefix = "data/angsd_pi/{ref}--{ssp}--{pop}"
    shell:
        """
        module load angsd
        angsd -GL 1 -P 5 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 30 -minQ 30 \
        -ref {input.ref}  -anc {input.trip} \
        -doMaf 2 -doMajorMinor 4 -doSaf 1 \
        -bam {input.bams} -out {params.prefix}
        """
rule pop_sfs:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        trip = "data/trip/trip_{ref}.fa.gz",
        bams = "data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt",
        saf = "data/angsd_pi/{ref}--{ssp}--{pop}.saf.gz"
    output:
        sfs = "data/angsd_pi/{ref}--{ssp}--{pop}.sfs"
    params:
        prefix = "data/angsd_pi/{ref}--{ssp}--{pop}",
    shell:
        """
        module load angsd
        realSFS {params.prefix}.saf.idx -P 10 > {params.prefix}.sfs
        realSFS saf2theta {params.prefix}.saf.idx -sfs {params.prefix}.sfs -outname {params.prefix}
        """

rule pop_pi:
    input:
        sfs = "data/angsd_pi/{ref}--{ssp}--{pop}.sfs"
    output:
        pi = "data/angsd_pi/{ref}--{ssp}--{pop}.{win}BP_theta.thetasWindow.gz.pestPG"
    params:
        prefix_in = "data/angsd_pi/{ref}--{ssp}--{pop}",
        prefix_out = "data/angsd_pi/{ref}--{ssp}--{pop}.{win}BP_theta",
        win = "{win}"
    shell:
        """
        module load angsd
        thetaStat do_stat {params.prefix_in}.thetas.idx -win {params.win} -step {params.win} -outnames {params.prefix_out}.thetasWindow.gz
        """

#thin
#awk 'NR ==1 || NR % 10 == 0' pi_test.txt  | wc -l  
#NGSadmix -likes angsd_til11_teo.beagle.gz -K 6 -P 10 -o ngsadmix_til11_teo -minMaf 0.05
