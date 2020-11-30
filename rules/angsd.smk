rule full_beagle:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        bams = "data/bamlist/{ref}-{ssp}-allBams.txt"
    output:
        "data/beagle/{ref}--{ssp}.beagle.gz"
    params:
        #prefix = "data/beagle/{ref}--{ssp}",
        scratch = my_scratch,
        prefix = my_scratch + "{ref}--{ssp}",
        final = "data/beagle/"
    shell:
        """
        mkdir -p {params.scratch}
    
        module load angsd 

        angsd -GL 1 -P 30 \
        -ref {input.ref} \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -minMapQ 30 -minQ 30 \
        -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam {input.bams} -out {params.prefix}

        mv {params.prefix}.beagle.gz {params.final}
        mv {params.prefix}.mafs.gz {params.final}
        mv {params.prefix}.arg {params.final}
 
        """

rule sample_gl:
    input:
        "data/beagle/{ref}--{ssp}.beagle.gz"
    output:
        "data/beagle/{ref}--{ssp}_100thin.beagle.gz"
    shell:
        """
        cat <(zcat {input} | head -n1) <(zcat {input} | tail -n+2 | sed -n '0~100p') | gzip > {output}
        """



rule admix:
    input:
        "data/beagle/{ref}--{ssp}_100thin.beagle.gz"
    output:
        "data/ngsAdmix/{ref}--{ssp}_K{K}.qopt"
    params:
        prefix = "data/ngsAdmix/{ref}--{ssp}_K{K}",
        K = "{K}"
    shell:
        """
        module load angsd
        NGSadmix -likes {input} -K {params.K} -P 10 -o {params.prefix} -minMaf 0.05
        """


#caution, this might fail because awk and snakemake are the pits
rule chico10_beagle:
    input:
        "data/beagle/{ref}--{ssp}_100thin.beagle.gz"
    output:
        "data/beagle/{ref}--{ssp}_100thin_random10_PalmarChico.beagle.gz"
    params:
        ssp = "{ssp}"
    shell:
        """
        zcat {input} | cut -f 1-3,`grep -Ff <(cat <(grep {params.ssp} pop_key | grep -v random | awk '{{print "Ind" NR-1 "\\t" $0}}' | cut -f1,4,5 | grep -v Palmar) <(grep {params.ssp} pop_key | grep -v random | awk '{{print "Ind" NR-1 "\\t" $0}}' | cut -f1,4,5 | grep Palmar | shuf -n 10 | sort -V) | awk '{{print $1"t"}}') <(zcat {input} | head -n1 | tr "\\t" "\\n" | cat -n | tail -n+4 | awk '{{print $1 "\\t" $2 "t"}}') | cut -f1 | tr "\\n" "," | sed 's/,$/\\n/g'` | gzip > {output}
        """

rule chico10_admix:
    input:
        "data/beagle/{ref}--{ssp}_100thin_random10_PalmarChico.beagle.gz"
    output:
        "data/ngsAdmix/{ref}_{ssp}_{k}_thin1M_random10_PalmarChico.qopt"
    params:
        prefix = "data/ngsAdmix/{ref}_{ssp}_{k}_thin1M_random10_PalmarChico",
        k = "{k}"
    shell:
        "NGSadmix -likes {input} -K {params.k} -P 10 -o {params.prefix} -minMaf 0.05"


rule chico_only_beagle:
    input:
        "data/beagle/{ref}--{ssp}_100thin.beagle.gz"
    output:
        "data/beagle/{ref}--{ssp}_100thin_PalmarChicoONLY.beagle.gz"
    params:
        ssp = "{ssp}"
    shell:
        """
        grep {params.ssp} pop_key | grep -v random | awk '{{print "Ind" NR-1 "\\t" $0}}' | grep  Palmar | cut -f1 > {params.ssp}PC_beagleID
        python src/subset_beagle.py -b {input} -i {params.ssp}PC_beagleID | gzip > {output}
        """
#"zcat {input} | cut -f1-3,124- | gzip > {output}"


rule chico_only_admix:
    input:
        "data/beagle/{ref}--{ssp}_100thin_PalmarChicoONLY.beagle.gz"
    output:
        "data/ngsAdmix/{ref}_{ssp}_{i}_thin1M_PalmarChicoONLY.qopt"
    params:
        prefix = "data/ngsAdmix/{ref}_{ssp}_{i}_thin1M_PalmarChicoONLY",
        k = "{i}"
    shell:
        """
        module load angsd
        NGSadmix -likes {input} -K {params.k} -P 10 -o {params.prefix} -minMaf 0.05
        """


rule split_beagle:
    input:
        "data/beagle/{ref}--{ssp}_100thin_PalmarChicoONLY.beagle.gz"
    output:
        "data/beagle/{ref}_{ssp}_{ch}_100thin_PalmarChicoONLY.beagle.gz"
    params:
        ref = "{ref}",
        ssp = "{ssp}",
        ch = "{ch}"
    shell:
        "cat <(zcat {input} | head -n1) <(zcat {input} | grep {params.ch}) | gzip > {output}"


rule chico_chrom:
    input:
        "data/beagle/{ref}_{ssp}_{ch}_100thin_PalmarChicoONLY.beagle.gz"
    output:
        "data/ngsAdmix/{ref}_{ssp}_{k}_{ch}_thin1M_PalmarChicoONLY.qopt"    
    params:
        ref = "{ref}",
        ssp = "{ssp}",
        prefix = "data/ngsAdmix/{ref}_{ssp}_{k}_{ch}_thin1M_PalmarChicoONLY",
        k = "{k}",
        ch = "{ch}"
    shell:
        """
        module load angsd
        NGSadmix -likes data/beagle/{params.ref}_{params.ssp}_{params.ch}_100thin_PalmarChicoONLY.beagle.gz -K {params.k} -P 10 -o {params.prefix} -minMaf 0.05
        """
        
        

#data/angsd_pi/v5--Teo--El_Rodeo.arg
#data/angsd_pi/v5--Teo--El_Rodeo.mafs.gz
#data/angsd_pi/v5--Teo--El_Rodeo.saf.gz
#data/angsd_pi/v5--Teo--El_Rodeo.saf.idx
#data/angsd_pi/v5--Teo--El_Rodeo.saf.pos.gz
#data/angsd_pi/v5--Teo--El_Rodeo.sfs    
#trip = "data/trip/trip_{ref}.fa.gz",

rule pop_beagle:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
        trip = "data/anc/{ref}_anc.fa",
        bams = "data/bamlist/{ref}--{ssp}--{pop}__bamlist.txt"
    output:
        saf = "data/angsd_pi/{ref}--{ssp}--{pop}.saf.gz"
    params:
        scratch = my_scratch,
        final = "data/angsd_pi/",
        prefix = my_scratch + "{ref}--{ssp}--{pop}"
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
        -ref {input.ref}  -anc {input.trip} \
        -doMaf 2 -doMajorMinor 4 -doSaf 1 \
        -bam {input.bams} -out {params.prefix}
        
         mv {params.prefix}.arg {params.final}
         mv {params.prefix}.mafs.gz {params.final}
         mv {params.prefix}.saf.gz {params.final}
         mv {params.prefix}.saf.idx {params.final}
         mv {params.prefix}.saf.pos.gz {params.final}
 
        """

rule pop_sfs:
    input:
        ref = "data/refs/{ref}/{ref}.fa",
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
