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

    std::cout << "Loading Model ... \n";
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

    std::cout << "Signatures Loaded \n";
}

void Signatures::predict(seqan::StringSet<seqan::Dna5String> &seqs, seqan::StringSet<seqan::CharString> &ids, std::vector<std::string> &readLabels, std::string &buffer, std::unordered_map<std::string, std::tuple<std::string, float>> &FuncPred)
{
    // TODO: This section has to be updated in the multiprocessing section'
    // Reads up to 10 records

    std::vector<std::string> readSeqs;
    std::vector<std::string> read_labels;
    std::vector<std::string> protSeqs;

    int num_reads = 1;
    int total_reads = 0;

    /////////////////////////////////////////////////////
    /*.............GET READING FRAMES.........*/
    /////////////////////////////////////////////////////

    // std::cout << "Open Reading Frames" << std::endl;
    // seqan::StringSet<seqan::String<seqan::AminoAcid>, seqan::Owner<seqan::ConcatDirect<>>> aaSeqs;
    // if (isreduced)
    // {
    //     seqan::GeneticCode<seqan::MURPHY10> GCode; // reduce the aminoacid alphabet to 10
    //     seqan::translate(aaSeqs, seqs, seqan::SIX_FRAME, GCode);
    // }
    // else
    // {
    //     seqan::GeneticCode<seqan::CANONICAL> GCode;
    //     seqan::translate(aaSeqs, seqs, seqan::SIX_FRAME, GCode);
    // }

    //////////////////////////////////////////////////////////////////////////////
    /*.............FILTER READS BY SIGNATURES AND CLASSIFY .........*/
    //////////////////////////////////////////////////////////////////////////////
    std::random_device rd;
    std::mt19937 rng(rd());
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
    // std::cout << "Traverse Reads file" << std::endl;

    // traverse the nucleotide sequences
    typedef seqan::Iterator<seqan::StringSet<seqan::String<seqan::AminoAcid>, seqan::Owner<seqan::ConcatDirect<>>>>::Type reading_frames_iterator;
    typedef seqan::Iterator<seqan::StringSet<seqan::DnaString>>::Type TStringSetIterator;

    seqan::StringSet<seqan::String<seqan::AminoAcid>, seqan::Owner<seqan::ConcatDirect<>>> reading_frames_sequences;
    seqan::GeneticCode<seqan::MURPHY10> GCode;

    for (unsigned dna_sequence = 0; dna_sequence < length(seqs); dna_sequence++)
    {
        // get the reading frames for each one
        seqan::translate(reading_frames_sequences, seqs[dna_sequence], seqan::SIX_FRAME, GCode);

        // traverse each reading frame
        for (reading_frames_iterator it = begin(reading_frames_sequences); it != end(reading_frames_sequences); ++it)
        {
            KMER.clear();

            // transform iterator to a string
            seqan::move(kimer, *it);
            toCSkmer = seqan::toCString(kimer);

            // length of the orf sequence
            l = toCSkmer.length();

            // select a random position
            std::uniform_int_distribution<int> uni(0, l - args->kmer - 1);
            rx = uni(rng);
            ishash = 0;
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

            // If there is at least one key in the hash table
            if (ishash > 0)
            {
                pre_buffer = KMER;

                // traverse the ORF and get the kmers for prediction
                // Using a sliding window of 1
                for (int ki = 0; ki < l - args->kmer - 1; ki++)
                {
                    pkmer = toCSkmer.substr(ki, args->kmer);
                    pre_buffer += ' ' + pkmer;
                    pkmer.clear();
                }

                buffer += pre_buffer + '\n';
                pre_buffer.clear();

                read_labels.push_back(seqan::toCString(ids[dna_sequence]));

                // store the sequences if they need to be reported in a file
                if (args->fastaOutput)
                {
                    std::stringstream iseq;
                    iseq << seqs[dna_sequence];
                    readSeqs.push_back(iseq.str());
                    iseq.clear();
                }
                else
                {
                    readSeqs.push_back(" ");
                }
            }
        }
    }

    std::stringstream trex(buffer);
    fasttext.predict(trex, 5, false, read_labels, 0, FuncPred, readSeqs, args->seq);

    // clearing variables to free memory
    trex.str(std::string());
    read_labels.clear();
    readSeqs.clear();
    buffer.clear();
    seqan::clear(seqs);
    seqan::clear(ids);
    seqan::clear(reading_frames_sequences);
}

void Signatures::Display(std::string message)
{
    std::cout << message << "\n";
}
