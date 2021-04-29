cat data/domestication_regions/domestication_regions.txt | tr "," "\t" | awk 'NR > 1 {print $1 "\t" $2-1 "\t" $3}' > data/domestication_regions/domestication_regions_v4.bed 

CrossMap.py bed data/chain_files/B73_RefGen_v4_to_Zm-B73-REFERENCE-NAM-5.0.chain data/domestication_regions/domestication_regions_v4.bed  | cut -f5- |  awk '{print "chr" $1 "\t" $2-25000 "\t" $3+25000}'  | bedtools merge -i stdin > data/domestication_regions/domestication_regions_v5.bed 


 paste data/domestication_regions/domestication_regions_v5.bed <(awk 'NR >1' data/domestication_regions/domestication_regions.txt) > data/domestication_regions/domestication_regions_v5_annotated.bed
