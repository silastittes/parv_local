StringVal="hard linkedHard linkedSoft neut soft"

cnn_pred () { for i in $StringVal; do python src/diploSHIC/diploSHIC.py predict --simData data/diploshic/cnns/$1.json  data/diploshic/cnns/$1.weights.hdf5 data/diploshic/$1/${i}.fvec ${i}_${1}.out; done; } 

pred_mat () { for i in $StringVal; do echo $i; tail -n+2 ${i}_${1}.out | cut -f1 | sort | uniq -c; done > ${1}_val.txt; } 

pred_emp() { python src/diploSHIC/diploSHIC.py predict data/diploshic/cnns/${1}.json  data/diploshic/cnns/${1}.weights.hdf5 $2 ${1}_empirical.preds; }

pred_table () { cut -f5 ${1}_empirical.preds | awk 'NR > 1' | sort | uniq -c  > ${1}_pred_summary.txt; }

run_preds() { 
    cnn_pred $1
    pred_mat $1
    pred_emp $1 $2
    pred_table ${1}

    for i in $StringVal; do rm ${i}_${1}.out; done
}

#run_preds v5--LR--Amatlan_de_Canas--chr1--0--308452471--s0.0001_0.0005--tau0.0_10000.0 data/diploshic/fvec_vcf/v5--LR--Amatlan_de_Canas--chr1--0--308452471.fvec 



