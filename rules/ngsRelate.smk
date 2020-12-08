rule ngsrelate_beagle:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        sites="data/refs/{ref}/{ref}_500M_sites.txt",
        sites_bin="data/refs/{ref}/{ref}_500M_sites.txt.bin",
        trip = "data/anc/{ref}_anc.fa",
        bams = "data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt"
    output:
        glf = "data/ngsRelate/{ref}--{ssp}--{pop}.glf.gz",
        maf = "data/ngsRelate/{ref}--{ssp}--{pop}.mafs.gz"
    params:
        scratch = my_scratch,
        prefix = my_scratch + "{ref}--{ssp}--{pop}",
        final = "data/ngsRelate/"
    shell:
        """
        mkdir -p {params.scratch}
        
        rm -f {params.prefix}.arg 
        rm -f {params.prefix}.glf.pos.gz
        rm -f {params.prefix}.mafs.gz
        rm -f {params.prefix}.glf.gz 

        module load angsd
        angsd -P 5 \
        -ref {input.ref} -anc {input.trip} \
        -sites {input.sites} \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minMapQ 30 -minQ 30 \
        -gl 1 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 \
        -b {input.bams} -out {params.prefix}

        mv {params.prefix}.arg {params.final} 
        mv {params.prefix}.glf.pos.gz {params.final}
        mv {params.prefix}.mafs.gz {params.final}
        mv {params.prefix}.glf.gz {params.final}                
        """

rule pop_freq:
    input:
        "data/ngsRelate/{ref}--{ssp}--{pop}.mafs.gz"
    output:
        "data/ngsRelate/{ref}--{ssp}--{pop}.freqs.txt"
    shell:
        "zcat {input} | cut -f7 | sed 1d > {output}"


rule ngsRelate:
    input:
        freq = "data/ngsRelate/{ref}--{ssp}--{pop}.freqs.txt",
        glf = "data/ngsRelate/{ref}--{ssp}--{pop}.glf.gz",
        bams = "data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt"
    output:
        "data/ngsRelate/{ref}--{ssp}--{pop}.ngsRelate.txt"
    shell:
        "src/ngsRelate/ngsRelate -p 5 -g {input.glf} -n `wc -l {input.bams} | cut -d ' ' -f1`  -f {input.freq}  -O {output}"
