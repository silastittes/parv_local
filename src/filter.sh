
vcftools --keep LR_key.txt --gzvcf data/vcf/v5/raw/v5.vcf.gz --max-alleles 2 --min-alleles 2 --remove-indels --minDP 5 --maxDP 100 --minGQ 30 --minQ 30 --max-missing 0.7 --recode --stdout | gzip > data/vcf/v5/filtered/v5_LR_filtered.vcf.gz &

vcftools --keep LR_key.txt --gzvcf data/vcf/til11/raw/til11.vcf.gz --max-alleles 2 --min-alleles 2 --remove-indels --minDP 5 --maxDP 100 --minGQ 30 --minQ 30 --max-missing 0.7 --recode --stdout | gzip > data/vcf/til11/filtered/til11_LR_filtered.vcf.gz &

vcftools --keep Teo_key.txt --gzvcf data/vcf/til11/raw/til11.vcf.gz --max-alleles 2 --min-alleles 2 --remove-indels --minDP 5 --maxDP 100 --minGQ 30 --minQ 30 --max-missing 0.7 --recode --stdout | gzip > data/vcf/til11/filtered/til11_Teo_filtered.vcf.gz &
        
vcftools --keep Teo_key.txt --gzvcf data/vcf/v5/raw/v5.vcf.gz --max-alleles 2 --min-alleles 2 --remove-indels --minDP 5 --maxDP 100 --minGQ 30 --minQ 30 --max-missing 0.7 --recode --stdout | gzip > data/vcf/v5/filtered/v5_Teo_filtered.vcf.gz 
        
bcftools query -l data/vcf/v5/filtered/v5_Teo_filtered.vcf.gz > data/vcf/v5/filtered/v5_Teo_filtered.key

bcftools query -l data/vcf/v5/filtered/v5_LR_filtered.vcf.gz > data/vcf/v5/filtered/v5_LR_filtered.key 

bcftools query -l data/vcf/til11/filtered/til11_Teo_filtered.vcf.gz > data/vcf/til11/filtered/til11_Teo_filtered.key

bcftools query -l data/vcf/til11/filtered/til11_LR_filtered.vcf.gz > data/vcf/til11/filtered/til11_LR_filtered.key


