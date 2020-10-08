grep Teo pop_key | awk '{print "data/bams/til11-alignments/deduped-bam/" $3 ".deduped.bam"}' > data/bamlist/til11-Teo-allBams.txt
grep LR pop_key | awk '{print "data/bams/til11-alignments/deduped-bam/" $3 ".deduped.bam"}' > data/bamlist/til11-LR-allBams.txt
grep Teo pop_key | awk '{print "data/bams/v5-alignments/deduped-bam/" $3 ".deduped.bam"}' > data/bamlist/v5-Teo-allBams.txt
grep LR pop_key | awk '{print "data/bams/v5-alignments/deduped-bam/" $3 ".deduped.bam"}' > data/bamlist/v5-LR-allBams.txt
