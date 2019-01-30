
#ifndef Signatures_H
#define Signatures_H

#include <iostream>
#include <random>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <thread>

#include "sparsepp/spp.h"
#include "hop/hopscotch_map.h"
#include <fasttext.h>
#include <sys/stat.h>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/translation.h>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/stream.h>

using spp::sparse_hash_map;

class Signatures
{

  public:
    std::string file_model;
    std::string fsignatures;

    fasttext::FastText fasttext;
    fasttext::FastText skipgram;

    std::vector<std::vector<std::string>> LABELS;
    std::vector<std::stringstream> KMER_LIST;

    tsl::hopscotch_map<std::string, bool> master_signature_hash;      // 7-mer
    tsl::hopscotch_map<std::string, bool> master_signature_hash_full; //11-mer

    int kmer_size;
    int seed_size;
    bool isreduced;

    std::shared_ptr<fasttext::Args> args;

    Signatures(std::shared_ptr<fasttext::Args>);
    void predict(seqan::StringSet<seqan::Dna5String> &seqs, seqan::StringSet<seqan::CharString> &ids, std::vector<std::string> &readLabels, std::string &buffer, std::unordered_map<std::string, std::tuple<std::string, float>> &FuncPred);
    void filter(seqan::StringSet<seqan::Dna5String> &seqs, seqan::StringSet<seqan::CharString> &ids, std::vector<std::string> &readLabels, std::string &buffer, std::unordered_map<std::string, std::tuple<std::string, float>> &FuncPred);
    void Display(std::string message);

  private:
    std::string MapSignatures(std::string header, std::string read);
};

#endif