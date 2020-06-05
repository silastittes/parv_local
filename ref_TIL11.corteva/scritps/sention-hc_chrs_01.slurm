#!/bin/bash
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -p freetest
#SBATCH -t 7-00:00:00
#SBATCH -J sention
#SBATCH -o sention.o%j
#SBATCH -e sention.e%j
#SBATCH --mail-user=arnstrm@gmail.com
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
IFS=$'\n' read -d '' -r -a CMDS < chrs_01
CMD=${CMDS[$SLURM_ARRAY_TASK_ID]}
echo "running sention on ${CMD}"
./runHaplotyper.sh ${CMD}
if [ $? -eq 0 ]
then
  echo "Success: sention on ${CMD}"
  exit 0
else
  echo "Failed: sention on ${CMD}"
exit 1
fi
