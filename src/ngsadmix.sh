module load angsd
genome=/group/jrigrp10/stittes/parv_local_v5_til11/data/refs/til11/Zea_mays_parviglumis_TIL11-inkind-chr.fasta
filters="-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0  -C 50  -minMapQ 20 -minQ 20"
bams="til11-Teo-allBams.txt"

#./angsd -GL 1 -out genolike -nThreads 10 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1  -bam bam.filelist

angsd -GL 1 -P 10 \
$filters \
-ref $genome -anc $genome \
-doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 \
-r chr1:0-500000000 \
-bam $bams -out angsd_til11_teo

NGSadmix -likes angsd_til11_teo.beagle.gz -K 6 -P 10 -o ngsadmix_til11_teo -minMaf 0.05
