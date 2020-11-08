
rule mop:
    input:
        cfile = "data/mop/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.txt",
        bams = "data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt"
    output:
        "data/mop/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.bed"
    params:
        chrom = "{chrom}",
        start = "{start}",
        end = "{end}"
    shell:
        "mop -m 0.7 -x 100 -i 5 -Q 30 -q 30 -b {input.bams} -R {params.chrom}:{params.start}-{params.end} > {output}"

rule mop_merge:
    input:
        #expand("data/mop/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.bed", zip, ref = mREF, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
        expand("data/mop/v5--{ssp}--{pop}--{chrom}--{start}--{end}.bed", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
    output:
        "data/mop/{ref}--{ssp}--{pop}_all.bed"
    run:
        shell("cat {input} > {output}")
