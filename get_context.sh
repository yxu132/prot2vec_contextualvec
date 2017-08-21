#!/bin/bash

#prot2vec_path=OUTPUT/uniprot-all_s100_w25_lr0.005_itr400_ns5_hs0/
query_sequences_fasta=$1
prot2vec_path=$2
window_size=$3
site_type=$4

python src/get_contexts.py $query_sequences_fasta $prot2vec_path --win_size=$window_size --site_type=$site_type