#include <iostream>
#include "Signatures.cpp"
#include "Index.cpp"
#include "args.h"
#include "utils.h"
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

//Standard namespace declaration
using namespace std;

void printUsage()
{
    std::cerr
        << "usage: metaMLP <command> <args>\n\n"
        << "The commands supported by metaMLP are:\n\n"
        << "  index              Index the reference sequences\n"
        << "  quant              Run quantification algorithm [prediction]\n"
        << std::endl;
}

void printPredictUsage()
{
    std::cerr
        << "usage: metaMLP quant <args>\n\n"
        << "The commands supported by metaMLP quant are:\n\n"
        << "  -model            Trained model - generated by metaMLP index\n"
        << "  -input            FASTA file with short sequence reads\n"
        << "  -output           output file to write the processed reads\n"
        << "  -kmer             k-mer size [default 11] same used during index\n\n"

        << "Optional Parameters: \n\n"

        << "  -tries            sensitivity [default 2] higher gives more hits more errors\n"
        << "  -minProbability   minimum probability to report sequences [default 0.8]\n"
        << "  -threads          number of threads to use\n"
        << "  -minReadChunkSize Load reads in memory [default 10000]\n\n"

        << "Optional Alphabet Parameters: \n\n"
        << "  -NoReduced        Enable it if index is built with the -NoReduced option\n"
        << std::endl;
}

void printIndexUsage()
{
    std::cerr
        << "usage: metaMLP index <args>\n\n"
        << "Mandatory Parameters:\n\n"
        << "  -input        Protein reference database\n"
        << "  -output       Output index\n"
        << "  -kmer         k-mer size (aminoacids) [default 11]\n"
        << "  -labp         Label index position in the FASTA header (default 4: >xx|xx|xx|label|xx)\n\n"

        << "Optional Training Parameters: \n\n"

        << "  -dim          word vector size [default 64] adjust for large number of classes\n"
        << "  -epoch        number of epochs for training the model [default 100]\n"
        << "  -lr           learning rate [default 1]\n"
        << "  -minn         minimum length of kmer ngram [default 3]\n"
        << "  -maxn         maximum length of kmer ngram [default 7]\n"
        << "  -ws           size of context window for building embeddings [default 5]\n"
        << "  -minCount     minimum ngram count [default 1]\n"
        << "  -wordNgrams   embedding joint ngrams [default 1]\n"
        << "  -loss         loss function {ns, hs, softmax} [default softmax]\n\n"

        << "Optional Alphabet Parameters: \n\n"
        << "  -NoReduced    Dissable the reduced alphabet and use all 20 Amino Acids\n\n"

        << "Optional Output: \n\n"
        << "  -fastaOutput    save fasta file with predicted sequences [default False]\n"

        << std::endl;
}

struct thread_data
{

    Signatures *signatures;
    std::string fi;
    // std::string report;
    int tid; // thread ID
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;
    std::vector<std::string> readLabels;
    std::string buffer;
    std::unordered_map<std::string, std::tuple<std::string, float>> FuncPred;
};

void thread_process(void *args)
{
    struct thread_data *params;
    params = (struct thread_data *)args;
    (*params->signatures).predict(params->seqs, params->ids, params->readLabels, params->buffer, params->FuncPred);
}

void compute_absolute_abundance()
{
}

void quant(int argc, char **argv)
{

    if (argc < 3)
    {
        printPredictUsage();
        exit(0);
    }

    std::shared_ptr<fasttext::Args> a = std::make_shared<fasttext::Args>();
    a->parseArgs(argc, argv);

    // std::cout << a->mink << std::endl;

    std::string report_file = a->output;

    int NUM_THREADS = a->thread;

    std::thread threads[NUM_THREADS];
    struct thread_data td[NUM_THREADS];

    // ************************************** //
    // how to use: ./map /path/to/fasta/file.fa /path/to/signatures.txt /path/to/fastx/model
    // int kmer_size = stoi(kmer);
    // Load the signatures in json format kmer-size
    Signatures signatures(a);

    // Load Fasta File
    std::ifstream input(a->input);
    seqan::SeqFileIn seqFileIn(seqan::toCString(a->input));
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

    // Open file and create Record Reader
    typedef seqan::Iterator<seqan::StringSet<seqan::Dna5String>>::Type NIter;
    seqan::CharString id;
    seqan::Dna5String seq;

    int ith = 0;
    int iseq = 0;
    int entries = 0;

    // ********************************************************************************************************
    // MAP SECTION
    // ********************************************************************************************************

    // read chunks of 10 reads at the time
    int chunks = a->minReadChunkSize;
    tsl::hopscotch_map<std::string, int> absolute_abundance;
    std::vector<std::string> label_sequence;
    int arglike = 0;

    std::ofstream fo(report_file);
    std::ofstream fabn(report_file + ".abn");
    std::ofstream fasta_o(report_file + ".fasta");
    if (!a->fastaOutput)
        remove(report_file.c_str());

    std::cout << "Processing Input File" << std::endl;

    while (!atEnd(seqFileIn))
    {
        for (int thread_number = 0; thread_number < NUM_THREADS; thread_number++)
        {
            seqan::readRecords(ids, seqs, seqFileIn, chunks);
            entries += length(ids);
            // move data to thread and clear
            seqan::move(td[thread_number].seqs, seqs);
            seqan::move(td[thread_number].ids, ids);
            td[thread_number].signatures = &signatures;
            td[thread_number].tid = thread_number;

            seqan::clear(seqs);
            seqan::clear(ids);

            // start a the new thread
            threads[thread_number] = std::thread(thread_process, &td[thread_number]);
        }

        for (int thread_number = 0; thread_number < NUM_THREADS; thread_number++)
        {
            threads[thread_number].join();
        }

        // display # processed reads
        std::cout << "processed reads " << entries << "\r" << std::flush;
    }

    // ********************************************************************************************************
    // REDUCE SECTION
    // ********************************************************************************************************

    for (int i = 0; i < NUM_THREADS; i++)
    {
        for (const auto &arglabel : td[i].FuncPred)
        {
            // if a minimum probability of 0.5 report the sequence
            if (std::get<1>(arglabel.second) >= a->minProbability)
            {
                // report individual classification
                label_sequence = fasttext::utils::splitString(std::get<0>(arglabel.second), '\t');

                // Report: sequence_id --> predicted_label --> probability
                fo << arglabel.first << "\t" << label_sequence[0] << "\t" << std::get<1>(arglabel.second) << "\n";

                // Report fasta file if enabled
                if (a->fastaOutput)
                {
                    fasta_o << ">" << arglabel.first << "|" << label_sequence[0] << std::endl
                            << label_sequence[1] << std::endl;
                }

                absolute_abundance[label_sequence[0]] += 1;
                arglike++;
            }
        }
    }

    std::cout << entries << " Processed sequences " << std::endl;
    std::cout << arglike << " Annotated sequences " << std::endl;

    // printout results:
    std::cout << "Computing Relative Abundances" << std::endl;
    std::stringstream ARGc;
    std::string HMP;
    for (const auto &item : absolute_abundance)
    {
        ARGc << item.first;
        HMP = ARGc.str();
        fabn << HMP.replace(HMP.length() - 2, HMP.length(), "").replace(0, 9, "") << "\t" << std::to_string(item.second) << std::endl;
        ARGc.str(std::string());
    }

    fabn.close();
    fo.close();

    exit(0);
}

void index(int argc, char **argv)
{
    if (argc < 3)
    {
        printIndexUsage();
        exit(0);
    }

    std::shared_ptr<fasttext::Args> a = std::make_shared<fasttext::Args>();
    a->parseArgs(argc, argv);
    Index index;
    index.indexing(a->input, a->output, a->kmer, a->labp - 1, a->reduced);

    // training model
    std::cout << "Indexing reference database ..." << std::endl;
    fasttext::FastText fastText;

    a->model = fasttext::model_name::sup;
    a->loss = fasttext::loss_name::hs;
    a->input = a->output + ".tr";
    // a->epoch = 100;
    // a->lr = 1;
    // a->minCount = 1;
    // a->tries = 5;
    // a->dim = 100;
    // a->wordNgrams = 2;
    fastText.train(a);

    std::cout << "Cleaning temporal files ..." << std::endl;
    // remove(a->input.c_str());

    exit(0);
}

//Main Function
int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printUsage();
        exit(0);
    }

    std::string command(argv[1]);
    if (command == "index")
    {
        index(argc, argv);
    }
    else if (command == "quant")
    {
        quant(argc, argv);
    }
    else
    {
        printUsage();
        exit(0);
    }

    return 0;
}
