/**
 * Copyright (c) 2016-present, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

#ifndef FASTTEXT_ARGS_H
#define FASTTEXT_ARGS_H

#include <istream>
#include <ostream>
#include <string>
#include "sparsepp/spp.h"
#include "hop/hopscotch_map.h"

namespace fasttext
{

enum class model_name : int
{
    cbow = 1,
    sg,
    sup
};
enum class loss_name : int
{
    hs = 1,
    ns,
    softmax
};

class Args
{
  public:
    Args();
    std::string input;
    std::string test;
    std::string output;
    bool reportSequences;
    bool fastaOutput;
    float minProbability;

    int minReadChunkSize;
    int kmer;
    int proc;
    int tries;
    int minSeqLen;
    std::string smodel;
    std::string db;
    int labp;
    int seed;
    int mink;
    bool reduced;

    double lr;
    int lrUpdateRate;
    int dim;
    int ws;
    int epoch;
    int minCount;
    int minCountLabel;
    int neg;
    int wordNgrams;
    loss_name loss;
    model_name model;
    int bucket;
    int minn;
    int maxn;
    int thread;

    double t;
    std::string label;
    int verbose;
    std::string pretrainedVectors;
    int saveOutput;

    bool qout;
    bool seq;
    bool retrain;
    bool qnorm;
    size_t cutoff;
    size_t dsub;

    // hash table for prediction
    tsl::hopscotch_map<std::string, bool> master_signature_hash_full;

    void parseArgs(int, char **);
    void printHelp();
    void save(std::ostream &);
    void load(std::istream &);
};

} // namespace fasttext

#endif
