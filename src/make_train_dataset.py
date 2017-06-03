# Load fasta file 

from Bio import SeqIO
import h5py
from operator import itemgetter
import json
from sklearn.feature_extraction.text import TfidfVectorizer
import sys

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
    }

# signatures = {i.split()[0]:i.split()[0] for i in open("../data/signatures_k11.txt")}

# USE:
# INPUT KMER-SIZE

k_mer_size = int(sys.argv[2]);
fo = open('../data/train_dataset.txt', 'w');
nop=0
print "processing input file .."

signatures = {}

for record in SeqIO.parse(open(sys.argv[1]),'fasta'):
    # Transform sequence to reduced aminoacids

    seq = ''.join([aa_index[ii] for ii in str(record.seq)])
    # print len(seq)
    # seq = str(record.seq);
    key = "__label__"+"".join(record.id.split('|')[3])+"__ "

    for s in range(0, len(seq)-k_mer_size,2):
        kml = []
        for ix, i in enumerate(range(s, len(seq)-k_mer_size,k_mer_size)):
            fragment =  seq[i:i+k_mer_size];
            signatures[fragment] = True;
            kml.append(fragment)
            if (ix+1)%5 == 0 :
                fo.write(key+" ".join(kml)+"\n")
                kml=[]

# print nop, 'genes without signatures'
fo.close();

print 'writing ',len(signatures),' signatures - kmers ...'
fo = open('../data/signatures.txt', 'w');

for i in signatures:
    fo.write(i+"\n");

fo.close()

print "training model ..."
cmd = "../bin/fasttext supervised -input ../data/train_dataset.txt -output ../data/AMR.model -epoch 50 -lr 1.0 -dim 10 -wordNgrams 2"
import os
os.system(cmd);

