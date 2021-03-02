params.slim_file = "$baseDir/dfe_alpha.slim"
params.prior_df = "$baseDir/prior_df.csv"

params_ch = Channel
    .fromPath(params.prior_df)
    .splitCsv(header: true)
    .map{ row -> [row.seed, row.sample_sizes, row.mu, row.c, row.loci, row.neg_mean, row.neg_shape, row.pos_mean, row.pos_shape, row.neg_prop, row.pos_prop, row.na, row.nb, row.n0, row.tb, row.t0]}
    //.view { row -> "${row.one} - ${row.two} - ${row.three}" }

process slim {
    
    input:
    tuple val(seed), val(sample_sizes), val(mu), val(c), val(loci), val(neg_mean), val(neg_shape), val(pos_mean), val(pos_shape), val(neg_prop), val(pos_prop), val(na), val(nb), val(n0), val(tb), val(t0) from params_ch

    output:
    file "${seed}_slim.txt" into sims_ch 

    """
    slim -s $seed \
    -define mu=$mu \
    -define c=$c \
    -define loci=$loci \
    -define neg_mean=$neg_mean \
    -define neg_shape=$neg_shape \
    -define pos_mean=$pos_mean \
    -define pos_shape=$pos_shape \
    -define neg_prop=$neg_prop \
    -define pos_prop=$pos_prop \
    -define na=$na \
    -define nb=$nb \
    -define n0=$n0 \
    -define tb=$tb \
    -define t0=$t0 \
    -define sample_sizes=\\"${sample_sizes}\\" $params.slim_file | grep line | cut -f2- > ${seed}_slim.txt 
    """
}

process combine {

    publishDir "$baseDir/simout/"

    input:
    file sims from sims_ch.collect()

    output:
    file "allsims.txt" into fullsims_ch


    """
    cat $sims > allsims.txt
    """

}

