#ifndef Index_H
#define Index_H

#include <string>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <fasttext.h>

#include <iostream>

using namespace std;

class Index{
    public:
        Index();
        void indexing(std::string finput, std::string output, int kmer, int label_index, bool isreduced);
        void training(std::string);
        std::string AA2Reduced(std::string aa);

    private:
        int valor;

};

#endif