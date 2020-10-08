module load angsd
genome=/group/jrigrp10/stittes/parv_local_v5_til11/data/refs/til11/Zea_mays_parviglumis_TIL11-inkind-chr.fasta
filters="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 20 -minQ 20"
bams="data/bamlist/til11--Teo--Palmar_Chico__bamlist.txt"

angsd -GL 1 -P 10 \
$filters \
-ref $genome  -anc $genome \
-doMaf 2 -doMajorMinor 4 -doSaf 1 \
-r chr1:10000000-13000000 \
-bam $bams -out chr1_10_20
realSFS chr1_10_20.saf.idx -P 10 -fold 1 > chr1_10_20.sfs
realSFS saf2theta chr1_10_20.saf.idx -sfs chr1_10_20.sfs -outname chr1_10_20
thetaStat do_stat chr1_10_20.thetas.idx -win 10000 -step 10000  -outnames chr1_10_20.theta.thetasWindow.gz
