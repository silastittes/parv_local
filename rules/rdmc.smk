
rule make_freq:
    input:
        "data/angsd_pi/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}.mafs.gz"
    output:
        "data/rdmc/freq/{ref}--{ssp}--{pop}--{chrom}--{start}--{end}_freq.bed.gz"
    shell:
        """
        zcat {input} | awk 'NR > 1 {{print $1 "\\t" $2-1 "\\t" $2 "\\t" $7}}' | gzip > {output}
        """

#FILTERED FOR MEAN ALLELE FREQUENCY BETWEEN 0.05 AND 0.95, WHICH IS A FAST AND HOPEFULLY A GOOD PROXY FOR VARIANCE
rule merge_freq:
    input:
        expand("data/rdmc/freq/{{ref}}--{pops}--{{chrom}}--{{start}}--{{end}}_freq.bed.gz", pops = FREQ_POPS)
    output:
        "data/rdmc/freq/{ref}--allpops--{chrom}--{start}--{end}_freq.txt.gz"
    shell:
        """
        bedtools unionbedg -i {input} -filler NA |\
        awk '{{for(i=4;i<=NF;i++) t+=$i; t_mean = t/(NF-3); if(t_mean> 0.05 && t_mean<0.95) print $0; t=0}}' | gzip > {output}
        """

#rule get_neutral:
#    input:
#        gff = "",
#        gbed = ""
#    output:
#    shell:
#zcat data/refs/v5/v5.gff3.gz | awk '$3 ~ /^gene$/ && $1 ~ /^chr/'  | bedtools slop -i stdin  -b 100000 -g data/refs/v5/v5.gbed | bedtools merge -i stdin | sort -V | bedtools complement -i stdin -g data/refs/v5/v5.gbed