import random
import sys


def writeToFile(list, file):
    print('Total number of sequences: '+str(len(list)))
    with open(file, 'w') as out:
        for line in list:
            out.write(line+'\n')

def gram_list(seq, gram):
    ret = []
    seq_size = len(seq)
    for i in range(0, seq_size, gram):
        if i + gram >= seq_size:
            ret.append(seq[i:])
        else:
            ret.append(seq[i:i+gram])
    return ret

def overlapping_gram_list(seq, gram):
    ret = []
    seq_size = len(seq)
    for i in range(0, seq_size-gram):
        ret.append(seq[i: i+gram])
    return ret

def uniprot_id_seq_map(file, name):
    prot_id = ''
    prot_seq = ''
    out = open('DATA/'+name+'_raw.txt', 'w')
    out_ids = open('DATA/'+name+'_id.txt', 'w')
    for line in open(file, 'r'):
        if line[0] == '>':
            if prot_id != '':
                out.write(prot_id + '\t' + prot_seq+'\n')
                out_ids.write(prot_id+'\n')
            prot_id = line.strip().split('|')[1]
            prot_seq = ''
        elif line.strip() != '':
            prot_seq = prot_seq + line.strip()

    if prot_id != '':
        out.write(prot_id + '\t' + prot_seq + '\n')
        out_ids.write(prot_id + '\n')


def protein_overlapping_3gram(name):
    lines = []
    with open('DATA/'+name+'_raw.txt', 'r') as input:
        lines.extend(input.readlines())
    examples = []
    for line in lines:
        seq = line.strip().split('\t')[1]
        examples.append(' '.join(overlapping_gram_list(seq, 3)))
    writeToFile(examples, 'DATA/'+name+'_overlap_3gram.txt')


def main(source_fasta_file):
    file_name = source_fasta_file.split('/')[-1]
    uniprot_id_seq_map(source_fasta_file, file_name[:-6])
    protein_overlapping_3gram(file_name[:-6])


if __name__ == '__main__':
    source_fasta_file = '/Volumes/YING_backup/uniprot-all.fasta'
    main(source_fasta_file)