rule maf2count:
    input:
        "data/angsd_pi/{ref}--{ssp}--{pop}.mafs.gz"
    output:
        "data/angsd_treemix/{ref}--{ssp}--{pop}.mafscount.bed"
    shell:
        """
        zcat < {input} | awk 'NR>1{{print $1 "\\t" $2-1 "\\t" $2 "\\t" int($7*$8) "," int((1-$7)*$8)}}' > {output}
        """

rule counts2treemix:
    input:
        expand("data/angsd_treemix/v5--{ssppop}.mafscount.bed", zip,  ssppop = mix_POP)
    output:
        "data/angsd_treemix/v5_treemix.gz"
    shell:
        "cat <(echo {input}) <(unionBedGraphs -i {input} | sed 's/\t0\t/\t0,0\t/g' | sed 's/\t0$/\t0,0/g' |  sed 's/\t0\t/\t0,0\t/g' | sed 's/\t0$/\t0,0/g' | cut -f4- ) | gzip > {output}"


rule filter_treemix:
    input:
        "data/angsd_treemix/{ref}_treemix.gz"
    output:
        "data/angsd_treemix/{ref}_treemix_filtered.gz"
    shell:
        "python src/filter_treemix.py -f {input} -m 0.05 -c 0.1 -t 100 | gzip > {output}"


rule treemix:
    input:
        "data/angsd_treemix/{ref}_treemix_filtered.gz"
    output:
        "data/angsd_treemix/{ref}_treemix.treeout.gz"
    params:
        prefix = "{ref}"
    shell:
        "src/treemix/src/treemix -i {input} -o {params.prefix}_treemix"
