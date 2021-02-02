rule make_freq:
    input:
        "data/angsd_pi/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.mafs.gz"
    output:
        "data/rdmc/freq/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}_freq.bed.gz"
    shell:
        """
        zcat {input} | awk 'NR > 1 {{print $1 "\\t" $2-1 "\\t" $2 "\\t" $7}}' | gzip > {output}
        """

#rule make_freq:
#    input:
#        "data/angsd_vcf/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.vcf.mop.gz"
#    output:
#        "data/rdmc/freq/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}_freq.bed.gz"
#    shell:
#        """
#        zcat {input} | grep -v "#"  | sed 's/;/\t/g' | awk '{{print $1 "\\t" $2-1 "\\t" $2 "\\t" $11}}'  | sed 's/AF=//g' | gzip > {output}
 #       """

#FILTERED FOR loci with zero variance
rule merge_freq:
    input:
        #expand("data/rdmc/freq/{{ref}}--{pops}--{{chrom}}--{{start}}--{{end}}_freq.bed.gz", pops = FREQ_POPS)
        expand("data/rdmc/freq/{pops}--{{chrom}}--{{start}}--{{end}}_freq.bed.gz", pops = FREQ_POPS)
    output:
        "data/rdmc/freq/{ref}--allpops--{chrom}--{start}--{end}_freq.txt.gz"
    shell:
        """
        bedtools unionbedg -i {input} -filler NA |\
        awk '{{for(i=4;i<=NF;i++) t+=$i; if(t > 0 && t < 1) print $0; t=0}}' |\
        grep -v NA | gzip > {output}
        """

#awk '{{for(i=4;i<=NF;i++) t+=$i; t_mean = t/(NF-3); if(t_mean> 0.05 && t_mean<0.95) print $0; t=0}}' |\
