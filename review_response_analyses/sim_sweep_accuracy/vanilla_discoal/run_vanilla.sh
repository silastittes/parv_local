size=50000


do_discoal(){
  ~/discoal/discoal 20 20 $size -t 500 -r 160 -wd 0 -a 500 -x 0.5 > discoal_out.txt
}

do_split(){
  
}

../../src/raisd-master/RAiSD -n TEST -I discoal_out.txt -L $size -f

csplit --digits=2  --quiet --prefix=outfile RAiSD_Report.TEST "/\/\//" "{*}"

do_bed(){
  for i in `ls outfile*`
  do
    awk '$2 > 5' $i | awk 'NR>1{print "ch1\t" $1-1 "\t" $1 "\t" $2}' | bedtools merge -i stdin -d 100000 -c 4 -o max > ${i}.bed
  done
  cd ..
}

do_bed
