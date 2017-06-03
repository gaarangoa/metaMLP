# Load fasta file 

from Bio import SeqIO
import h5py
from operator import itemgetter
import json
from sklearn.feature_extraction.text import TfidfVectorizer
import sys

k_mer_size = int(sys.argv[2]);

aa_index = {
    'Y': 'F',
    'F': 'F',
    'L': 'I',
    'I': 'I',
    'V': 'I',
    'M': 'I',
    'C': 'C',
    'W': 'F',
    'D': 'E',
    'N': 'E',
    'T': 'S',
    'S': 'S',
    'Q': 'E',
    'K': 'K',
    'E': 'E',
    'R': 'K',
    'A': 'A',
    'G': 'G',
    'H': 'H',
    'P': 'P',
    'X': 'X',
    'Z': 'X',
};

print "loading input file ..."

db = {}
for record in SeqIO.parse(open(sys.argv[1]),'fasta'):
    # print record
    seq = "".join([aa_index[ik] for ik in str(record.seq)]);
    key = record.id.split('|')[3]
    kml = []
    for i in range(0, len(seq)-k_mer_size,1):
        kml.append(seq[i:i+k_mer_size])
    try:
        db[key].append(" ".join(kml))
    except:
        db[key] = [" ".join(kml)]

# TF - IDF
ptop = 1;
print 'Computing TF-IDF, where the top ',ptop*100,'% of kmers are selected ... '

signatures = {}

for type in db:
    D = db[type]
    # 
    tfidf = TfidfVectorizer(use_idf=True,
                                norm="l2", 
                                min_df = 1,
                                ngram_range=(1,1),
                                sublinear_tf=False)
    # 
    A = tfidf.fit_transform(D);
    # 
    idf = tfidf._tfidf.idf_
    IDF = zip(tfidf.get_feature_names(), idf)
    # 
    top = int(len(IDF)*ptop);
    signatures[type] = [i[0].upper() for i in sorted(IDF, key=itemgetter(1))[:top]]

S = {}

for i in signatures:
    for j in signatures[i]:
        S[j]='';

print 'Storing signatures ...'

fo = open('../data/signatures_k11.txt', 'w')
for i in S:
    # fo.write("".join([aa_index[k] for k in i])+"\n");
    fo.write("".join([k for k in i])+"\n");

# json.dump(S, open('../map_signatures/signatures.k-5AA.json','w'));

# json.dump(signatures, open('../signatures/signatures.k-5AA.json','w'))



