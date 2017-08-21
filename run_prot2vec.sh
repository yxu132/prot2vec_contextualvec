#!/bin/bash

source_database_path=$1

echo "Step 1: Pre-processing: splitting sequences into 3-gram biological words..."
python src/input_data.py $source_database_path

echo "Step 2: Traning: prot2vec models based on the given protein sequence database... (It may take days. )"
python src/train_prot2vec.py $source_database_path _overlap_3gram --size=100 --window=25 --train_alpha=0.005 --iter=400 --negative=5 --hs=0 --infer_epoch=30000 --infer_alpha=0.01

