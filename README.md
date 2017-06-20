## PredX
PredX is a fast approach to annotate antimicrobial resistance from metagenomics short sequence reads. It compares sequence reads against a protein database, similar to what blastx does. However, instead of performing sequence alignment, it uses a linear classifier using word embeddings (kmer embeddings) to represent kmers as vectors. We use the FastText library to perform text classification and the seqan library to perform the six reading frame translation. To index protein sequences (reference database) we use the reduced alphabet proposed by RAPSearch. PredX is several orders of magnituf fast than blastx, diamond and deepARG. However, it is not designed to annotate low-similarity sequences (deepARG wins here!!).

### Status
In development-testing close to release, and you shouldn't be able to see this repository!!

## Installation
ARGfast does not require any external library, so, it can be installed by just typing: 
        
        git clone https://github.com/gaarangoa/ARGfast.git
        cd ARGfast
        make

## Quick start guide
predX can be used to annotate other functional categories apart of Antibiotic Resistance. To do so, you need to follow the next steps to index a reference database. 

### Full command line options
#### predX
        usage: predX <command> <args>

        The commands supported by predX are:

        index              Index the reference sequences
        quant              Run quantification algorithm [fasttext model]

#### Index
        usage: predX index <args>

        The commands supported by predX index are:

        -input        Protein reference database
        -output       Output index
        -kmer         k-mer size in aminoacids [default 11]
        -labp         Label index position in the FASTA header (default 4: >xx|xx|xx|label|xx)
        -NoReduced    Dissable the reduced alphabet and use all 20 Amino Acids

#### Annotation
        usage: predX quant <args>

        The commands supported by predX quant are:

        -model        Trained model - generated by predX index
        -input        FASTA file with short sequence reads
        -output       output file to write the processed reads
        -proc         Number of threads to use [default 8]
        -kmer         k-mer size [default 11 - same kmer used in predX index]
        -seed         seed size [default 11 - amino acids]
        -mink         minimum number of kmers that each read has to contain [default 5]
        -NoReduced    Enable it if index is built with the -NoReduced option

To create the index of the Antibiotic Resistance profiling, you can use our developed database [1] (please check the publication). 

### Index protein reference file
To index a protein reference database you first need to add the gene labels of the database to the header of the fasta file. For instance, just add the label to one of the fields in the header separated by the character "|" as follows >gene_id|field1|label|field3. For this particular example the index where the label is placed is 3. 

        ./predX index -input /path/to/reference_protein_database.fasta -output /path/to/output/data/prefix -labp 3

Where -labp corresponds to the index where the label is placed within the fasta header and prefix is the database index name. predX generates two files prefix.bin and prefix.kh. 


### Perform Functional Profiling
To perform the functional profiling you need to use the quant command. The input file has to be a fasta file, it does not support fastq files, to get a fasta file you can use another tool such as vsearch to merge paired end reads. The -model depicts the model created by ./predX index.

        ./predX quant -input /path/to/input/reads.fasta  -model /path/to/prefix -output /path/to/out/report.txt

## Antibiotic Resistance
The predX program was designed specially to perform a fast antimicrobial resistance profiling. You can download the precomputed antibiotic resistance index from here and jump to the ./predX quant command. 