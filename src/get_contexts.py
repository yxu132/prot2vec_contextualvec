import numpy as np
import gensim.models as model
import os, sys


class Config():
    infer_epoch = 30000
    infer_alpha = 0.01
    window_size = 3
    site_type = '*'
    def __init__(self, args):
        for par in args:
            parts = par.split('=')
            if '--infer_epoch' in par:
                Config.infer_epoch = int(parts[1])
            elif '--infer_alpha' in par:
                Config.infer_alpha = float(parts[1])
            elif '--win_size' in par:
                Config.window_size = int(parts[1])
            elif '--site_type' in par:
                Config.site_type = parts[1]


def readLine(path):
    ret = []
    with open(path, 'r') as file:
        ret.extend(file.readlines())
    ret = [line.strip() for line in ret]
    return ret

def read_windows(aas, window, filter=True, type='Y'):
    print 'TYPE:\t', type
    ret = []
    remains = []
    for aa in aas:
        remain = []
        vec = []
        for ind in xrange(len(aa)):
            if filter:
                if type == '*':
                    if aa[ind] == 'S' or aa[ind] == 'Y' or aa[ind] == 'T':
                        start = max([0, ind-window])
                        end = min([len(aa), ind+window+1])
                        vec.append(aa[start:end])
                        remain.append(ind)
                else:
                    if aa[ind] == type:
                        start = max([0, ind-window])
                        end = min([len(aa), ind+window+1])
                        vec.append(aa[start:end])
                        remain.append(ind)
            else:
                start = max([0, ind - window])
                end = min([len(aa), ind + window + 1])
                vec.append(aa[start:end])
                remain.append(ind)
        ret.append(vec)
        remains.append(remain)
    return np.array(ret), remains

def input_str(input):
    ids, seqs = [], []
    input = input.strip()
    if '\n' in input:
        lines = input.split('\n')
        current_id, current_seq = '', ''
        for line in lines:
            if line.startswith('>'):
                if current_id != '':
                    ids.append(current_id)
                    seqs.append(current_seq)
                current_id = line.strip()[1:]
                current_id = current_id.split('|')[1]
                current_seq = ''
            else:
                current_seq += line.strip()
        if current_id != '':
            ids.append(current_id)
            seqs.append(current_seq)
        else:
            return ids, seqs
    else:
        if ', ' in input:
            ids = input.split(', ')
        elif ',' in input:
            ids = input.split(',')
        else:
            ids.append(input)
        for i in xrange(len(ids)):
            seqs.append('')
    return ids, seqs

class Prot2vecModel:

    def __init__(self, input_model_path, infer_alpha, infer_epoch):

        # print input_model_path
        self._model = model.Doc2Vec.load(input_model_path)
        self.infer_alpha = infer_alpha
        self.infer_epoch = infer_epoch

    def ngram(self, seq):
        ret = []
        for i in range(0, len(seq) - 3):
            ret.append(seq[i: i + 3])
        return ret

    def infer(self, n_grams, alpha, steps):
        prot2vec = self._model.infer_vector(n_grams, alpha=alpha, steps=steps)
        return prot2vec

    def get(self, seq):
        n_grams = self.ngram(seq)
        prot2vec = self.infer(n_grams, alpha=self.infer_alpha, steps=self.infer_epoch)
        return prot2vec


def vecs_to_string(vecs):
    ret = ''
    for vec in vecs:
        item = '['
        for i in vec:
            item += str(i) + ' '
        ret += item + ']\n'
    return ret

def get(args):
    input = args[1]
    prot2vec_path = args[2]
    cfg = Config(args[3:])
    print 'INPUT:\t', input
    input = readLine(input)
    input = '\n'.join(input)
    ids, seqs = input_str(input)

    if len(ids) > 0:
        windows_, remains = read_windows(seqs, cfg.window_size, type=cfg.site_type)
    else:
        return 'Invad input. '

    if os.path.exists('CACHE/ids_'+str(cfg.window_size)+'.txt'):
        inferred_ids = readLine('CACHE/ids_'+str(cfg.window_size)+'.txt')
        inferred_ids_origin_len = len(inferred_ids)
        inferred_vecs = list(np.load('CACHE/vecs_'+str(cfg.window_size)+'.npy'))
    else:
        inferred_ids = []
        inferred_ids_origin_len = 0
        inferred_vecs = []

    map_ = dict(zip(inferred_ids, inferred_vecs))

    ptm_model = Prot2vecModel(prot2vec_path+'/overlap_3gram/prot2vec.bin', cfg.infer_alpha, cfg.infer_epoch)

    ret = 'position \t aa \t context \t generated vector \n'
    for ind, wins in enumerate(windows_):
        ret += '>'+ids[ind]+'\n'
        remain = remains[ind]
        for indj, win in enumerate(wins):
            if win in map_:
                prot2vec = map_[win]
            else:
                prot2vec = ptm_model.get(win)
                map_[win] = prot2vec
                inferred_ids.append(win)
                inferred_vecs.append(prot2vec)
            prot2vec = [str(item) for item in list(prot2vec)]
            ret += str(remain[indj]+1) + '\t' + str(seqs[ind][remain[indj]]) + '\t' + win + '\t' + ','.join(prot2vec) + '\n'

    if len(inferred_ids) > inferred_ids_origin_len:
        with open('CACHE/ids_'+str(cfg.window_size)+'.txt', 'w') as output:
            output.write('\n'.join(inferred_ids))
        np.save('CACHE/vecs_'+str(cfg.window_size)+'.npy', np.array(inferred_vecs))

    return ret

def main(args):
    ret = get(args)
    file_name = args[1].split('/')[-1]
    file_name = file_name[:file_name.index('.')]
    with open(file_name+'.out', 'w') as output:
        output.write(ret)
    print 'Execution completed. '

if __name__ == '__main__':
    main(sys.argv)


