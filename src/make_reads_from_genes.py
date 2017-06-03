import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def chunkstring(string, length):
    pieces = range(0, len(string), length)+range(length/2, len(string), length) + range(length/3, len(string), length) + range(length/4, len(string), length) + range(length/5, len(string), length);
    return (string[0+i:length+i] for i in pieces)

fi=sys.argv[1] # fastafile from besthit
fo=sys.argv[2]
tag=sys.argv[3]
read_size = int(sys.argv[4])

reads=[]

for record in SeqIO.parse(fi, "fasta"):
    seqn = str(record.seq);
    id = record.id.replace('|FEATURES', '');
    l=0;
    chunks = chunkstring(seqn, read_size);
    for k in chunks:
        if len(k)>=read_size-10:
            rseq = Seq(k.upper());
            rid = tag+"_"+id+'|'+str(l);
            reads.append(SeqRecord(rseq, id=rid, name='', description=''));
            l+=1;

SeqIO.write(reads, open(fo, 'w'), 'fasta')