#include "Signatures.h"
#include "functions.h"

using spp::sparse_hash_map;
using namespace std;

std::mutex mtx;

// TODO!!: this split function is repeated in the INDEX.cpp, I have to move it to some common library so i can used it multiple times.
template <typename Out>
void splitx(const std::string &s, char delim, Out result)
{
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        *(result++) = item;
    }
}

std::vector<std::string> splitx(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    splitx(s, delim, std::back_inserter(elems));
    return elems;
}
// END:TODO!!

// Constructor
// Load the signatures from the input file
Signatures::Signatures(std::shared_ptr<fasttext::Args> a)
{

    args = a;

    file_model = a->smodel;
    fsignatures = file_model + ".kh";
    kmer_size = a->kmer;
    seed_size = a->seed;
    isreduced = a->reduced;

    /////////////////////////////////////////////////////
    /*.............LOADING FASTTEXT MODEL.........*/
    /////////////////////////////////////////////////////

    std::cout << "Loading Fast Text Model ... \n";
    fasttext.loadModel(file_model + ".bin");

    //////////////////////////////////////////////////////
    /*.............LOADING SIGNATURES TO MEMORY.........*/
    //////////////////////////////////////////////////////

    std::cout << "Loading Signatures ... \n";

    std::ifstream ifs(fsignatures);
    std::string line;
    std::vector<std::string> iline;

    while (std::getline(ifs, line))
    {
        iline = splitx(line, '\t');
        master_signature_hash_full[iline[0]] = true; //iline[1];
    }
    ifs.close();
}

void Signatures::predict(seqan::StringSet<seqan::Dna5String> &seqs, seqan::StringSet<seqan::CharString> &ids, std::vector<std::string> &readLabels, std::string &buffer, std::unordered_map<std::string, std::tuple<std::string, float>> &FuncPred)
{
    // TODO: This section has to be updated in the multiprocessing section'
    // Reads up to 10 records

    std::vector<std::string> readSeqs;
    std::vector<std::string> protSeqs;

    int num_reads = 1;
    int total_reads = 0;

    /////////////////////////////////////////////////////
    /*.............GET READING FRAMES.........*/
    /////////////////////////////////////////////////////
    seqan::StringSet<seqan::String<seqan::AminoAcid>, seqan::Owner<seqan::ConcatDirect<>>> aaSeqs;
    if (isreduced)
    {
        seqan::GeneticCode<seqan::MURPHY10> GCode; // reduce the aminoacid alphabet to 10
        seqan::translate(aaSeqs, seqs, seqan::SIX_FRAME, GCode);
    }
    else
    {
        seqan::GeneticCode<seqan::CANONICAL> GCode;
        seqan::translate(aaSeqs, seqs, seqan::SIX_FRAME, GCode);
    }

    typedef seqan::Iterator<seqan::StringSet<seqan::String<seqan::AminoAcid>, seqan::Owner<seqan::ConcatDirect<>>>>::Type AIter;

    //////////////////////////////////////////////////////////////////////////////
    /*.............FILTER READS BY SIGNATURES AND CLASSIFY .........*/
    //////////////////////////////////////////////////////////////////////////////
    std::random_device rd;
    std::mt19937 rng(rd());
    // std::uniform_int_distribution<int> uni(0, int(read_length/3)-seed_size-1);
    int frame = 1, rx;
    typedef seqan::Infix<seqan::String<seqan::AminoAcid>>::Type kmer;

    seqan::String<char> kimer;
    std::string toCSkmer;
    std::string KMER;
    std::string pre_buffer;

    std::size_t stop_c;

    std::string pkmer;
    std::string pseed;

    // possible kmer sizes
    std::uniform_int_distribution<int> random_mer_size(3, args->kmer);
    int variable_kmer = random_mer_size(rng);

    int l = 0;
    int ishash = 0;
    // mtx.lock();
    int iframe = 0;
    for (AIter it = begin(aaSeqs); it != end(aaSeqs); ++it)
    {

        KMER.clear();

        if (frame == 6)
        {
            frame = 1;
            total_reads++;
        }
        else
        {
            frame++;
        }

        std::uniform_int_distribution<int> uni(0, length(*it) - args->kmer - 2);
        rx = uni(rng); // random position
        // Move iterator to a string-like structure

        seqan::move(kimer, *it);
        toCSkmer = seqan::toCString(kimer);

        // If the read has a stop codon, go to next reading frame:
        stop_c = toCSkmer.find_first_of('*');
        if (stop_c < 30)
            continue;

        l = toCSkmer.length();

        // Get kmer from a random position

        ishash = 0;
        // make n tries to get the right kmer from the read
        int tries = 0;
        while (1)
        {
            rx = uni(rng);
            KMER = toCSkmer.substr(rx, args->kmer);
            ishash = master_signature_hash_full.count(KMER);
            if (ishash > 0)
                break;
            if (tries == args->tries)
                break;
            tries++;
        }

        int manykmers = 0;
        // Got a kmer at all?, great, make a sentence and predict!! :)
        if (ishash > 0)
        {
            pre_buffer = KMER;

            for (int ki = 0; ki < l - args->kmer; ki += args->kmer)
            {
                pkmer = toCSkmer.substr(ki, args->kmer);
                pre_buffer += ' ' + pkmer;
                pkmer.clear();
            }

            if (length(ids[total_reads]) > 1)
            {
                // Check if the read has a proper header
                buffer += pre_buffer + '\n';
                pre_buffer.clear();

                readLabels.push_back(seqan::toCString(ids[total_reads]));
                std::stringstream iseq;
                iseq << seqs[total_reads];
                if (args->seq)
                {
                    readSeqs.push_back(iseq.str());
                }

                num_reads++;
            }
        }

        iframe++;
    }

    seqan::clear(seqs);
    seqan::clear(aaSeqs);
    seqan::clear(ids);

    std::stringstream trex(buffer);
    buffer.clear();
    fasttext.predict(trex, 1, false, readLabels, 0, FuncPred, readSeqs, args->seq);

    trex.str(std::string());
    readLabels.clear();
    readSeqs.clear();

    seqan::clear(aaSeqs);
}

void Signatures::Display(std::string message)
{
    std::cout << message << "\n";
}
