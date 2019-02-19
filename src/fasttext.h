/**
 * Copyright (c) 2016-present, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree. An additional grant
 * of patent rights can be found in the PATENTS file in the same directory.
 */

#ifndef FASTTEXT_FASTTEXT_H
#define FASTTEXT_FASTTEXT_H

#define FASTTEXT_VERSION 11 /* Version 1a */
#define FASTTEXT_FILEFORMAT_MAGIC_INT32 793712314

#include <time.h>

#include <atomic>
#include <memory>
#include <set>

#include "args.h"
#include "dictionary.h"
#include "matrix.h"
#include "qmatrix.h"
#include "model.h"
#include "real.h"
#include "utils.h"
#include "vector.h"
// #include "Signatures.h"

#include <sparsepp/spp.h>
#include <hop/hopscotch_map.h>
#include <sys/stat.h>

using spp::sparse_hash_map;

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/translation.h>
#include <seqan/find.h>
#include <seqan/index.h>
#include <seqan/stream.h>
namespace fasttext
{

class FastText
{
  private:
  public:
    FastText();

    std::shared_ptr<Args> args_;
    std::shared_ptr<Dictionary> dict_;

    std::shared_ptr<Matrix> input_;
    std::shared_ptr<Matrix> output_;

    std::shared_ptr<QMatrix> qinput_;
    std::shared_ptr<QMatrix> qoutput_;

    std::shared_ptr<Model> model_;

    // std::shared_ptr<Signatures> signatures_;

    std::atomic<int64_t> tokenCount;

    clock_t start;
    void signModel(std::ostream &);
    bool checkModel(std::istream &);

    bool quant_;

    sparse_hash_map<std::string, int> absolute_abundance;

    std::string fout;
    typedef seqan::Iterator<seqan::StringSet<seqan::String<seqan::AminoAcid>, seqan::Owner<seqan::ConcatDirect<>>>>::Type reading_frames_iterator;

    void setOutput(std::string fout);
    void setReadLabel(std::string label, int tid);
    void Report(std::string fo);

    void getVector(Vector &, const std::string &);
    void saveVectors();
    void saveOutput();
    void saveModel();
    void loadModel(std::istream &);
    void loadModel(const std::string &);
    void printInfo(real, real);

    void supervised(Model &, real, const std::vector<int32_t> &,
                    const std::vector<int32_t> &);
    void cbow(Model &, real, const std::vector<int32_t> &);
    void skipgram(Model &, real, const std::vector<int32_t> &);
    std::vector<int32_t> selectEmbeddings(int32_t) const;
    void quantize(std::shared_ptr<Args>);
    void test(std::istream &, int32_t);
    void predict2(std::istream &, int32_t, bool, std::vector<std::string>, int, std::unordered_map<std::string, std::tuple<std::string, float>> &, std::vector<std::string>, bool);
    void predict_line(
        std::istream &,
        int32_t,
        std::string,
        std::ofstream &);

    void predict(
        std::istream &,
        int32_t,
        std::vector<std::pair<real, std::string>> &) const;

    void wordVectors();
    void sentenceVectors();
    void ngramVectors(std::string);
    void textVectors();
    void printWordVectors();
    void printSentenceVectors();
    void precomputeWordVectors(Matrix &);
    void findNN(const Matrix &, const Vector &, int32_t,
                const std::set<std::string> &);
    void nn(int32_t, std::string);
    void analogies(int32_t);
    void trainThread(int32_t);
    void train(std::shared_ptr<Args>);

    void mapThread(int32_t);
    void map(std::shared_ptr<Args>);

    void loadVectors(std::string);

    void process_read(
        std::string header,
        std::string sequence,
        std::ofstream &ofile,
        seqan::Dna5String dna_sequence,
        seqan::GeneticCode<seqan::MURPHY10> GCode,
        seqan::StringSet<seqan::String<seqan::AminoAcid>, seqan::Owner<seqan::ConcatDirect<>>> reading_frames_sequences);
};

} // namespace fasttext
#endif
