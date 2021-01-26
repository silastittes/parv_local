#!/bin/bash

snakemake --jobs 50 \
    --rerun-incomplete \
    --latency-wait 60 \
    --cluster-config submit.json  \
    --cluster "sbatch -p {cluster.p} --mem {cluster.mem} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --job-name {cluster.name} -e {cluster.e} -o {cluster.o}"

snakemake --jobs 50 \
    --snakefile RDMC.smk \
    --rerun-incomplete \
    --latency-wait 60 \
    --cluster-config submit.json  \
    --cluster "sbatch -p {cluster.p} --mem {cluster.mem} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --job-name {cluster.name} -e {cluster.e} -o {cluster.o}"


