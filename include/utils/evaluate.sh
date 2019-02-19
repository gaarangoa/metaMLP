

grep ">" test_reads.fasta | awk '{i+=1; split($0, x, "|"); print i"\t"x[3];}' > test_reads.labels
paste ../data/test_reads.fasta.ARG.index ../data/test_reads.fasta.ARG.abn | awk '{gsub("__label__","",$0); gsub("__","",$0); print}' > performance.txt

python ../src/utils/evaluate.py ../data/test_reads.labels ../data/performance.txt 


