import gensim.models as model
import gensim.matutils as matutils
from gensim.models.doc2vec import TaggedDocument
import numpy as np
import sys, logging, os

class Config():
    size = 100
    window = 25
    min_count = 1
    sample = 1e-5
    workers = 16
    hs = 0
    dm = 0
    negative = 5
    dbow_words = 1
    dm_concat = 0
    iter = 300
    infer_epoch = 20000
    infer_alpha = 0.01
    train_alpha = 0.025
    train_min_alpha = 0.0001
    def __init__(self, args):
        for par in args:
            parts = par.split('=')
            if '--size' in par:
                Config.size = int(parts[1])
            elif '--window' in par:
                Config.window = int(parts[1])
            elif '--min_count' in par:
                Config.min_count = int(parts[1])
            elif '--sample' in par:
                Config.sample = float(parts[1])
            elif '--negative' in par:
                Config.negative = int(parts[1])
            elif '--dbow_words' in par:
                Config.dbow_words = int(parts[1])
            elif '--iter' in par:
                Config.iter = int(parts[1])
            elif '--infer_epoch' in par:
                Config.infer_epoch = int(parts[1])
            elif '--infer_alpha' in par:
                Config.infer_alpha = float(parts[1])
            elif '--train_alpha' in par:
                Config.train_alpha = float(parts[1])
            elif '--train_min_alpha' in par:
                Config.train_min_alpha = float(parts[1])
            elif '--hs' in par:
                Config.hs = int(parts[1])

def infer(input_doc, topN, input_model, output, cfg):
    documents = []
    for ind, line in enumerate(open(input_doc, 'r').readlines()):
        if ind < topN:
            documents.append(line.strip().split())
    m = model.Doc2Vec.load(input_model+'/prot2vec.bin')
    infer_res = []
    for ind, doc in enumerate(documents):
        print('Inferring document: '+str(ind))
        infer_res.append(m.infer_vector(doc, alpha=cfg.infer_alpha, steps=cfg.infer_epoch))
    np.save(output, infer_res)

def train_prot2vec(input, output, cfg):
    documents = []
    if 'interpro' in input:
        for line in open(input, 'r'):
            comps = line.strip().split('\t')
            print comps
            document = TaggedDocument(words=comps[1].split(' '), tags=[comps[0]])
            documents.append(document)
    else:
        documents = model.doc2vec.TaggedLineDocument(input)
    m = model.Doc2Vec(documents, size=cfg.size, window=cfg.window,
                      alpha=cfg.train_alpha, min_alpha=cfg.train_min_alpha,
                      min_count=cfg.min_count, sample=cfg.sample,
                      workers=cfg.workers, hs=cfg.hs, dm=cfg.dm,
                      negative=cfg.negative, dbow_words=cfg.dbow_words,
                      dm_concat=cfg.dm_concat, iter=cfg.iter)
    m.save(output+'/prot2vec.bin')

def eval_similarty_n(v1, v2):
    vv1 = [v for v in v1]
    vv2 = [v for v in v2]
    sims = []
    for ind, v_1 in enumerate(vv1):
        v_2 = vv2[ind]
        sim = np.dot(matutils.unitvec(np.array(v_1)), matutils.unitvec(np.array(v_2)))
        print('similarity '+str(ind)+': '+str(sim))
        sims.append(sim)
    print('avg similarity: '+str(np.mean(sims, axis=0)))

def main(args):
    file_name = args[1].split('/')[-1]
    train_path = 'DATA/'+ file_name[:-6] + args[2] + '.txt'

    # set up configuration parameters
    cfg = Config(args[3:])
    logging.basicConfig(format='%(asctime)s: %(levelname)s: %(message)s', level=logging.INFO)

    # output directory
    output_dir = 'OUTPUT/'+file_name[:-6]+'_s'+str(cfg.size)+'_w'+str(cfg.window)+'_lr'+str(cfg.train_alpha)+'_itr'+str(cfg.iter)+'_ns'+str(cfg.negative)+'_hs'+str(cfg.hs)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    doc2vec_dir = output_dir + '/' + args[2][1:]
    if not os.path.exists(doc2vec_dir):
        os.mkdir(doc2vec_dir)

    # traning
    train_prot2vec(train_path, doc2vec_dir, cfg)


if __name__ == "__main__":
    # command = 'python train_prot2vec.py ../DATA/UNIPROT _non_overlap_3gram --size=100 --window=25 --train_alpha=0.005 --iter=400 --negative=5 --hs=0 --infer_epoch=30000 --infer_alpha=0.01'
    main(sys.argv)
