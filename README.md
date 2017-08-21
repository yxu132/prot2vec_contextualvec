# prot2vec_contextualvec

This repository is for training prot2vec models and generating contextual vectors for S/T/Y sites in protein sequences. 

## Basic requirement

To run the scripts in this repository, you need to install Python >= 2.7.x, Numpy >= 1.3, SciPy >= 0.7, Gensim >= 2.3.0.  

## Training your own prot2vec

To train your own prot2vec model, you will need to prepare a database of protein sequences in the format of FASTA. Say the path to your protein sequence database is $path_to_database, you can train a prot2vec model by running the following script, 

> ./run_prot2vec.sh $PATH_TO_DATABASE 

This training process could take several days depending on the size of your protein sequence database. 

The trained prot2vec model will be under the directory DATA.  

## Obtain the contextual vectors for protein sequences

After the training of the prot2vec model is completed, you can obtain the contextual feature vectors for S/T/Y sites from your own prot2vec model. 

The protein sequences in query should be saved in a file in the format of FASTA. If the path to this query file is $PATH_TO_QUERY, you can obtain the contextual vectors (window size = 2*3+1) for Y sites by running the following script, 

> get_context.sh $PATH_TO_QUERY $PATH_TO_PROT2VEC 3 Y

Here, the $PATH_TO_PROT2VEC refers to the path where your prot2vec model is saved, you can find it in directory DATA. The window size can be adjusted to any size that you would like to test and the site type can be one of ('S', 'T', 'Y' and '*') where * means all three types of sites. 



