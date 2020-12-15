import os

mSSP, mPOP, mCHROM, mSTART, mEND = glob_wildcards("data/mop/v5--{ssp}--{pop}--{chrom}--{start}--{end}.txt")
mop_files = expand("data/mop/v5--{ssp}--{pop}--{chrom}--{start}--{end}.bed", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
chroms = list(set(mCHROM))

#configuration 
mu = 3e-8
diploshic="src/diploSHIC/diploSHIC.py"
discoal="src/discoal/"
sim_bps = 55000
site_idx = list(range(11))
sweep_locs = [0.045454545454545456, 0.13636363636363635, 0.22727272727272727, 0.3181818181818182, 0.4090909090909091, \
    0.5, 0.5909090909090909, 0.6818181818181818, 0.7727272727272727, 0.8636363636363636, 0.9545454545454546]
tau = [[1e4, 5e4], [0.0, 1e4]]
alpha = [[5e-5, 1e-4], [1e-4, 5e-4], [5e-4, 1e-3]]

def build_discoal(discoal_path, discoal_file, sweep_locs, sim_bps, tau, alpha, mu, out_path, out_prefix):
    with open(discoal_file) as d:
        discoal_string = discoal_path + d.readline().strip()

    #grab theta and N0
    theta_low, theta_high = discoal_string.split("-Pt ")[1].split(" ")[0:2]
    theta = (float(theta_low) + float(theta_high))/2
    N0 = theta / (4*mu*sim_bps)
    
    #write commands
    for i,x in enumerate(sweep_locs):
        for a in alpha:
            for t in tau:
                s1 = 2*N0*a[0]
                s2 = 2*N0*a[1]
                tmid = ((t[0] + t[1])/2)/(4*N0)
                t1 = t[0]/(4*N0)
                t2 = t[1]/(4*N0)
                neutral_string = f"{discoal_string} -x {x}"
                hard_string = f"{neutral_string} -ws {tmid} -Pa {s1} {s2} -Pu {t1} {t2}"
                soft_string = f"{hard_string} -Pf 0 0.1"
                file_name = f"{out_path}/TYPE--{out_prefix}--s{a[0]}_{a[1]}--tau{t[0]}_{t[1]}--window_{i}.sh"
                with open(file_name.replace("TYPE", "NEUTRAL"), "w") as n:
                    print(f"{neutral_string}", file = n)
                with open(file_name.replace("TYPE", "HARD"), "w") as h:
                    print(f"{hard_string}", file = h)
                with open(file_name.replace("TYPE", "SOFT"), "w") as s:
                    print(f"{soft_string}", file = s)
rule all:
    input:
        [f"data/diploshic/vcf_pred/v5--Teo--random1_Palmar_Chico--chr1--0--308452471--s{s[0]}_{s[1]}--tau{t[0]}_{t[1]}.pred" for s in alpha for t in tau],
        [f"data/diploshic/vcf_pred/v5--LR--random1_Palmar_Chico--chr1--0--308452471--s{s[0]}_{s[1]}--tau{t[0]}_{t[1]}.pred" for s in alpha for t in tau],
        [f"data/diploshic/fvec_vcf/v5--Teo--random2_Palmar_Chico--chr1--0--308452471--s{s[0]}_{s[1]}--tau{t[0]}_{t[1]}.fvec" for s in alpha for t in tau]
        #[f"data/diploshic/cnns/v5--Teo--random1_Palmar_Chico--chr1--0--308452471--s{s[0]}_{s[1]}--tau{t[0]}_{t[1]}.weights.hdf5" for s in alpha for t in tau]
        #[f"data/diploshic/discoal/SOFT--v5--Teo--random1_Palmar_Chico--chr1--0--308452471--s{s[0]}_{s[1]}--tau{t[0]}_{t[1]}--window_1.out.gz" for s in alpha for t in tau]
        #"data/diploshic/cnns/v5--Teo--random1_Palmar_Chico--chr1--0--308452471--s0_0.01--tau0_0.05.weights.hdf5"

rule build_discoal:
    input:
        command = "data/mushi/{ref}--{ssp}--{pop}--mushi_discoal.txt"
    output:
        expand("data/diploshic/discoal/{{type}}--{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}--window_{window}.sh", window = site_idx)
    params:
        prefix = "{ref}--{ssp}--{pop}--{chrom}--{start}--{end}",
        path = "data/diploshic/discoal/"
    run:
        build_discoal(discoal, input.command, sweep_locs, sim_bps, tau, alpha, mu, params.path, params.prefix)

rule discoal:
    input:
        "data/diploshic/discoal/{type}--{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}--window_{window}.sh"
    output:
        "data/diploshic/discoal/{type}--{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}--window_{window}.out.gz"
    run:
        shell("bash {input} | gzip > {output}")
    
rule mask:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        gbed = "data/refs/{ref}/{ref}.gbed",
        bed = "data/mop/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.bed"
    output:
        mask_chrom = temp("data/diploshic/mask/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--mask.fa"),
        fasta_chrom = temp("data/diploshic/mask/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--raw.fa")
    params:
        chrom = "{chrom}"
    shell:
        """
        samtools faidx {input.ref} {params.chrom} > {output.fasta_chrom}
        bedtools complement -g {input.gbed} -i {input.bed} | bedtools maskfasta -fi {output.fasta_chrom} -bed stdin -fo {output.mask_chrom} 
        """
#bedtools maskfasta -fi {output.fasta_chrom} -bed {input.bed} -fo {output.mask_chrom} 
rule fvecsim:
    input:
        mask = "data/diploshic/mask/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--mask.fa",
        discoal = "data/diploshic/discoal/{type}--{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}--window_{window}.out.gz"
    output:
        "data/diploshic/fvec/{type}--{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}--window_{window}.fvec"
    params:
        chrom = "{chrom}"
    shell:
        """
        python {diploshic} fvecSim diploid {input.discoal} {output} --totalPhysLen {sim_bps} --maskFileName {input.mask} --chrArmsForMasking {params.chrom}
        """


rule trainingsets:
    input:
        neut = expand("data/diploshic/fvec/NEUTRAL--{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}--window_{window}.fvec", zip, window = site_idx),
        neut_sim = "data/diploshic/fvec/NEUTRAL--{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}--window_0.fvec",
        soft = expand("data/diploshic/fvec/SOFT--{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}--window_{window}.fvec", zip, window = site_idx),
        hard = expand("data/diploshic/fvec/HARD--{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}--window_{window}.fvec", zip, window = site_idx),
    output:
        neut = "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/neut.fvec",
        soft = "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/soft.fvec",
        hard = "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/hard.fvec",
        linkedHard = "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/linkedHard.fvec",
        linkedSoft = "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/linkedSoft.fvec"
    params:
        simdir = "data/diploshic/fvec/",
        outdir =  "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/"
    shell:
        """
        mkdir -p {params.outdir}
        python {diploshic} makeTrainingSets {input.neut_sim} \
               {params.simdir}SOFT \
               {params.simdir}HARD \
                5 \
                0,1,2,3,4,6,7,8,9,10 \
               {params.outdir}
        """
#data/diploshic/cnns/v5--Teo--random1_Palmar_Chico--chr1--0--308452471--s{s[0]}_{s[1]}--tau{t[0]}_{t[1]}.weights.hdf5
rule train:
    input:
        rules.trainingsets.output
    output:
        wt = "data/diploshic/cnns/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}.weights.hdf5",
        json = "data/diploshic/cnns/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}.json",
    params:
        prefix = "data/diploshic/cnns/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}",
        outdir = "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/" 
    shell:
        "python {diploshic} train {params.outdir} {params.outdir} {params.prefix}"


rule fvec_vcf:
    input:
        vcf = "data/angsd_vcf/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.vcf.gz",
        mask = "data/diploshic/mask/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--mask.fa"
    output:
        "data/diploshic/fvec_vcf/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}.fvec"
    params:
        chrom = "{chrom}",
        end = "{end}"
    shell:
        "python src/diploSHIC/diploSHIC.py fvecVcf diploid {input.vcf} {params.chrom} {params.end} {output} --winSize {sim_bps}  --maskFileName {input.mask}"

rule emp_predict:
    input:
        wt = "data/diploshic/cnns/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}.weights.hdf5",
        json = "data/diploshic/cnns/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}.json",
        fvec = "data/diploshic/fvec_vcf/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}.fvec"
    output:
        "data/diploshic/vcf_pred/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}.pred"
    shell:
        "python src/diploSHIC/diploSHIC.py predict {input.json} {input.wt} {input.fvec} {output}"



include: "../rules/demography.smk"

