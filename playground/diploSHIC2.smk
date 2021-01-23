import pickle

mSSP, mPOP, mCHROM, mSTART, mEND = glob_wildcards("data/mop/v5--{ssp}--{pop}--{chrom}--{start}--{end}.txt")
mop_files = expand("data/mop/v5--{ssp}--{pop}--{chrom}--{start}--{end}.bed", zip, ssp = mSSP, pop = mPOP, chrom = mCHROM, start = mSTART, end = mEND)
chroms = list(set(mCHROM))

#configuration 
mu = 3e-8
c = 1.6e-8 
diploshic="src/diploSHIC/diploSHIC.py"
discoal="src/discoal/discoal"
sim_bps = 110000
n_sims = 2000
site_idx = list(range(11))
sweep_locs = [0.045454545454545456, 0.13636363636363635, 0.22727272727272727, 0.3181818181818182, 0.4090909090909091, \
    0.5, 0.5909090909090909, 0.6818181818181818, 0.7727272727272727, 0.8636363636363636, 0.9545454545454546]

#not really tau, just generations
tau = [[1e4, 5e4], [0.0, 1e4]]
#not really alpha, just s
alpha = [[5e-5, 1e-4], [1e-4, 5e-4], [5e-4, 1e-3]]


##testing!!!!!!!!
#not really tau, just generations
tau = [[0.0, 1e4]]
#not really alpha, just s
alpha = [[5e-4, 1e-3]]



#theta_0 = 4*N_0*mu*args.discoal_nsites
#rho = 4*N_0*c*args.discoal_nsites
#with open(f"{args.prefix}_discoal.txt", "w") as discoal_file:
#    print(f"discoal {ns} {args.discoal_reps} {args.discoal_nsites} -Pt {0.8*theta_0} {1.2*theta_0} -Pre {rho} {3*rho} {en_string}", file = discoal_file)

def build_discoal(discoal_path, discoal_file, sim_type, sweep_locs, sim_bps, tau, alpha, mu, out_path, out_prefix):
    with open(discoal_file, "rb") as d:
        discoal_dict = pickle.load(d)
        
    #grab theta and N0
    N_0 = discoal_dict['N_0']
    ns = discoal_dict['ns']
    en_string = discoal_dict['en_string']
    theta_0 = 4*N_0*mu*sim_bps
    rho = 4*N_0*c*sim_bps
    discoal_string = f"{discoal} {ns} {n_sims} {sim_bps} -Pt {0.8*theta_0} {1.2*theta_0} -Pre {rho} {3*rho} {en_string}"
    
    alpha_f = alpha.split("_")
    tau_f = tau.split("_")
    
    s1 = 2*N_0*float(alpha_f[0])
    s2 = 2*N_0*float(alpha_f[1])
    t1 = float(tau_f[0])/(4*N_0)
    t2 = float(tau_f[1])/(4*N_0)
    tmid = (t1 + t2)/2
  
    #write commands
    for i,x in enumerate(sweep_locs):
        neutral_string = f"{discoal_string} -x {x}"
        hard_string = f"{neutral_string} -ws {tmid} -Pa {s1} {s2} -Pu {t1} {t2}"
        soft_string = f"{hard_string} -Pf 0 0.1"
        file_name = f"{out_path}/{sim_type}--{out_prefix}--s{alpha}--tau{tau}--window_{i}.sh"

        if sim_type == "NEUTRAL":
            with open(file_name, "w") as s:
                print(f"{neutral_string}", file = s)

        if sim_type == "SOFT":
            with open(file_name, "w") as s:
                print(f"{soft_string}", file = s)

        if sim_type == "HARD":
            with open(file_name, "w") as s:
                print(f"{hard_string}", file = s)


#data/diploshic/vcf_pred/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.pred
rule all:
    input:
        [f"data/diploshic/vcf_pred/v5--Teo--random1_Palmar_Chico--chr1--0--308452471--s{s[0]}_{s[1]}--tau{t[0]}_{t[1]}/v5--Teo--random1_Palmar_Chico--chr1--0--308452471--{0}--{1000000}.pred" for s in alpha for t in tau]


rule build_discoal:
    input:
        command = "data/mushi/{ref}--{ssp}--{pop}--mushi_discoal.txt"
    output:
        expand("data/diploshic/discoal/{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}/{{type}}--{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}--window_{window}.sh", window = site_idx)
    params:
        sim_type = "{type}",
        tau = "{tau}",
        alpha = "{s}",
        prefix = "{ref}--{ssp}--{pop}--{chrom}--{start}--{end}",
        path = "data/diploshic/discoal/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/"
    run:
        build_discoal(discoal, input.command, params.sim_type, sweep_locs, sim_bps, params.tau, params.alpha, mu, params.path, params.prefix)

rule discoal:
    input:
        "data/diploshic/discoal/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/{type}--{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}--window_{window}.sh"
    output:
        "data/diploshic/discoal/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/{type}--{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}--window_{window}.out.gz"
    run:
        shell("bash {input} | gzip > {output}")
    
rule mask:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        gbed = "data/refs/{ref}/{ref}.gbed",
        bed = "data/mop/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.bed"
    output:
        mask_chrom = "data/diploshic/mask/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--mask.fa",
        fasta_chrom = "data/diploshic/mask/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--raw.fa"
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
        discoal = "data/diploshic/discoal/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/{type}--{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}--window_{window}.out.gz"
    output:
        "data/diploshic/fvec/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/{type}--{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}--window_{window}.fvec"
    params:
        stat_dir = "data/diploshic/fvec/{type}--{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}_stats",
        chrom = "{chrom}"
    shell:
        """
        python {diploshic} fvecSim diploid {input.discoal} {output} --outStatsDir {params.stat_dir} --totalPhysLen {sim_bps} --maskFileName {input.mask} --chrArmsForMasking {params.chrom}
        """

rule trainingsets:
    input:
        neut = expand("data/diploshic/fvec/{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}/NEUTRAL--{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}--window_{window}.fvec", zip, window = site_idx),
        neut_sim = "data/diploshic/fvec/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/NEUTRAL--{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}--window_0.fvec",
        soft = expand("data/diploshic/fvec/{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}/SOFT--{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}--window_{window}.fvec", zip, window = site_idx),
        hard = expand("data/diploshic/fvec/{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}/HARD--{{ref}}--{{ssp}}--{{pop}}--{{chrom}}--{{start}}--{{end}}--s{{s}}--tau{{tau}}--window_{window}.fvec", zip, window = site_idx),
    output:
        neut = "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/neut.fvec",
        soft = "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/soft.fvec",
        hard = "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/hard.fvec",
        linkedHard = "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/linkedHard.fvec",
        linkedSoft = "data/diploshic/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/linkedSoft.fvec"
    params:
        simdir = "data/diploshic/fvec/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/",
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
        fvec = "data/diploshic/fvec_vcf/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--{vstart}--{vend}.fvec",
        statfile = "data/diploshic/fvec_vcf/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--{vstart}--{vend}_stats.txt"
    params:
        chrom = "{chrom}",
        end = "{end}",
        vstart = "{vstart}",
        vend = "{vend}"
    shell:
        """
        python src/diploSHIC/diploSHIC.py fvecVcf diploid {input.vcf} {params.chrom} {params.end} {output.fvec} \
            --winSize {sim_bps} --maskFileName {input.mask} --segmentStart {params.vstart} --segmentEnd {params.vend} --statFileName {output.statfile}
        """

#in order to train, mask file from a specific chromosome is needed. In order NOT to train a separate net for every chromosome, tell snakemake which net to use on which chromosome.
#name output directory according the net, file according to which chromosome the *empirical* data came frome.
rule emp_predict:
    input:
        wt = "data/diploshic/cnns/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}.weights.hdf5",
        json = "data/diploshic/cnns/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}.json",
        fvec = "data/diploshic/fvec_vcf/{eref}--{essp}--{epop}--{echrom}--{estart}--{eend}--{vstart}--{vend}.fvec"
    output:
        "data/diploshic/vcf_pred/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}--s{s}--tau{tau}/{eref}--{essp}--{epop}--{echrom}--{estart}--{eend}--{vstart}--{vend}.pred"
    shell:
        "python src/diploSHIC/diploSHIC.py predict {input.json} {input.wt} {input.fvec} {output}"

include: "../rules/demography.smk"

