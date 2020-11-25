rule ngsrelate_beagle:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        #trip = "data/trip/trip_{ref}.fa.gz",
        trip = "data/anc/{ref}_anc.fa",
        bams = "data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt"
    output:
        glf = "data/ngsRelate/{ref}--{ssp}--{pop}.glf.gz",
        maf = "data/ngsRelate/{ref}--{ssp}--{pop}.mafs.gz"
    params:
        prefix = "data/ngsRelate/{ref}--{ssp}--{pop}"
    shell:
        """
        module load angsd
        angsd -P 5 \
        -ref {input.ref} -anc {input.trip} \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minMapQ 30 -minQ 30 \
        -gl 1 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 \
        -b {input.bams} -out {params.prefix}
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
