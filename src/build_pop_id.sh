run_pop(){
    num_id=$1
    jrial=$2
    pop_n=$3
    subspp=$4
    pop=$5
    for (( c=1; c<=$num_id; c++ )); do echo -e "$pop_n\t$subspp\t$jrial-${c}\t$5"; done
}
#EXAMPLE
run_pop 55 "JRIAL1" "pop_5" "LR" Palmar_Chico

run_pop 50 "JRIAL8" "pop_6" "Teo" Palmar_Chico
