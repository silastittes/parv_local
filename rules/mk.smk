rule fold_nuc_sites:
    input:
        maf = "data/angsd_pi/{ref}--{ssp}--{pop}.mafs.gz",
        fold = "data/trip/trip_v5_FOLD{fold}",
    output:
        "data/angsd_sfs/{ref}--{ssp}--{pop}_fold{fold}_{nucs}_sites.txt"
    params:
        nucs = "{nucs}"
    shell:
        """
        module load angsd
        python src/parse_fold_type_maf.py -m {input.maf} -f {input.fold} -n {params.nucs} > {output}
        """

rule site_index:
    input:
        "data/angsd_sfs/{ref}--{ssp}--{pop}_fold{fold}_{nucs}_sites.txt"
    output:
        "data/angsd_sfs/{ref}--{ssp}--{pop}_fold{fold}_{nucs}_sites.txt.bin"
    shell:
        """
        module load angsd
        angsd sites index {input}
        """


rule sfs_beagle:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        trip = "data/trip/trip_{ref}.fa.gz",
        bams = "data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt",
        sites_bin = "data/angsd_sfs/{ref}--{ssp}--{pop}_fold{fold}_{nucs}_sites.txt.bin",
        sites = "data/angsd_sfs/{ref}--{ssp}--{pop}_fold{fold}_{nucs}_sites.txt"
    output:
        saf = "data/angsd_sfs/{ref}--{ssp}--{pop}--fold{fold}_{nucs}.saf.idx"
    params:
        prefix = "data/angsd_sfs/{ref}--{ssp}--{pop}--fold{fold}_{nucs}"
    shell:
        """
        module load angsd
        angsd -GL 1 -P 15 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 30 -minQ 30 \
        -sites {input.sites} \
        -ref {input.ref}  -anc {input.trip} \
        -doMaf 2 -doMajorMinor 4 -doSaf 1 \
        -bam {input.bams} -out {params.prefix}
        """

rule fold_nuc_sfs:
    input:
        #saf = "data/angsd_pi/{ref}--{ssp}--{pop}.saf.idx",
        saf = "data/angsd_sfs/{ref}--{ssp}--{pop}--fold{fold}_{nucs}.saf.idx",
    output:
        "data/angsd_sfs/{ref}--{ssp}--{pop}_fold{fold}_{nucs}_sfs.txt"
    shell:
        """
        module load angsd
        realSFS -P 15 {input.saf} > {output}
        """
