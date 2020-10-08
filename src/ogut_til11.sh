src/hal/bin/halLiftover  --outPSL data/TIL11-B73_halfiles/Zea_mays_parviglumis-TIL11_B73.hal B73 <(awk '$1 ~ /[0-9]/' data/map/ogut_v5.map.txt | awk '{print "chr" $1 "\t" $2 "\t" $2 + 1 "\t" $3}') Zea_mays_parviglumis-TIL11 T

awk 'BEGIN{print "chr\tpos\tcm"}{print $1 "\t" $2 "\t" $4}' T > data/map/ogut_til11.map.txt 

rm T
