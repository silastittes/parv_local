# SNP calling for JRIAL lines

The sentieon based pipeline for calling variants.

## Input data

Data was downloaded from CyVerse.

```
/iplant/home/rossibarra/rare_alleles_seq/JRIAL11/JRIAL11_fastq
/iplant/home/rossibarra/rare_alleles_seq/LR/JRIAL1/JRIAL1_fastq
/iplant/home/rossibarra/rare_alleles_seq/JRIAL10/JRIAL10_fastq
/iplant/home/rossibarra/rare_alleles_seq/teosinte_50/JRIAL8/JRIAL8_fastq
```

> there was a problem with JRIAL8_10-13 lines and was re-done using the source files provided by JRI.


## Alignment and QC

Each read pairs was analysed with the [`runSention.sh`](scripts/runSention.sh) script. Simple files with list of reads ([`reads_00`](assets/reads_00) and [`reads_01`](assets/reads_01)) was created and was processed using slurm array jobs ([`sention-hc_chrs_00.slurm`](scritps/sention-hc_chrs_00.slurm) and [`sention-hc_chrs_01.slurm`](scritps/sention-hc_chrs_01.slurm)) .

```
sbatch --array=0-97 sention-hc_chrs_00.slurm
sbatch --array=0-96 sention-hc_chrs_01.slurm
```

Once the jobs completed, the reports and bam files were separated and uploaded to CyVerse.

```bash
mkdir -p deduped-bam sorted-bam reports jobfiles
mkdir reports/{LOG,METRICS,PDF,RECAL-DATA}
mv *deduped.bam* deduped-bam/
mv *sorted.bam* sorted-bam/
mv *.pdf reports/PDF/
mv *.csv *.post reports/RECAL-DATA
mv *_metrics.txt *summary.txt reports/METRICS
```

Upload:

```bash
imkdir /iplant/home/aseetharam/TIL11-alignments
icd /iplant/home/aseetharam/TIL11-alignments
iput -a -K -P -r deduped-bam/
iput -a -K -P -r reports/
```

## Variant calling

The `senteion` was used for running the `HaplotypeCaller`. The genome was split into chunks are were processed independently.


```bash
genome=
ml picard
ml samtools
ml bwa
ml bedtools2
ml bioawk
picard CreateSequenceDictionary REFERENCE=${genome} OUTPUT=${genome%.*}.dict
samtools faidx ${genome}
bwa index ${genome}
bioawk -c fastx '{print $name"\t"length($seq)}' ${genome} | grep "^chr" > ${genome%.*}.length
bedtools makewindows -w 16000000 -g ${genome%.*}.length | awk '{print $1":"$2+1"-"$3}' > ${genome%.*}.chunks.txt
split -d -l ${genome%.*}.chunks.txt chrs_
```

These genome chunks were processed using [`runHaplotyper.sh`](scripts/runHaplotyper.sh) script, using the slurm files [`sention-hc_chrs_00.slurm`](scritps/sention-hc_chrs_00.slurm) and [`sention-hc_chrs_01.slurm`](scritps/sention-hc_chrs_01.slurm).


The results (vcf files) were uploaded to CyVerse.
