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


