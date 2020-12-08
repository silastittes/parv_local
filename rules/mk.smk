rule fold_nuc_sites:
    input:
        maf = list(set(expand("data/angsd_pi/{{ref}}--{{ssp}}--{{pop}}--{chrom}--{start}--{end}.mafs.gz", zip, chrom = mCHROM, start = mSTART, end = mEND))),
        #fold = "data/trip/trip_v5_FOLD{fold}",
        #fold = "data/refs/{ref}/{ref}_FOLD{fold}",
        fold = "data/anc/{ref}_anc_FOLD{fold}"
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
        trip = "data/anc/{ref}_anc.fa",
        bams = "data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt",
        sites_bin = "data/angsd_sfs/{ref}--{ssp}--{pop}_fold{fold}_{nucs}_sites.txt.bin",
        sites = "data/angsd_sfs/{ref}--{ssp}--{pop}_fold{fold}_{nucs}_sites.txt"
    output:
        saf = "data/angsd_sfs/{ref}--{ssp}--{pop}--fold{fold}_{nucs}.saf.idx"
    params:
        scratch = my_scratch,
        prefix = my_scratch + "{ref}--{ssp}--{pop}--fold{fold}_{nucs}",
        final = "data/angsd_sfs/"
    shell:
        """
        mkdir -p {params.scratch}
        
        rm -f {params.prefix}.arg
        rm -f {params.prefix}.mafs.gz
        rm -f {params.prefix}.saf.gz
        rm -f {params.prefix}.saf.idx
        rm -f {params.prefix}.saf.pos.gz

        module load angsd
        angsd -GL 1 -P 15 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 30 -minQ 30 \
        -sites {input.sites} \
        -ref {input.ref}  -anc {input.trip} \
        -doMaf 2 -doMajorMinor 4 -doSaf 1 \
        -bam {input.bams} -out {params.prefix}

         mv {params.prefix}.arg {params.final}
         mv {params.prefix}.mafs.gz {params.final}
         mv {params.prefix}.saf.gz {params.final}
         mv {params.prefix}.saf.idx {params.final}
         mv {params.prefix}.saf.pos.gz {params.final}
        """

rule fold_nuc_sfs:
    input:
        saf = "data/angsd_sfs/{ref}--{ssp}--{pop}--fold{fold}_{nucs}.saf.idx",
    output:
        "data/angsd_sfs/{ref}--{ssp}--{pop}_fold{fold}_{nucs}_sfs.txt"
    shell:
        """
        module load angsd
        realSFS -P 15 {input.saf} > {output}
        """
