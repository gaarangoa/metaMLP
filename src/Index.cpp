#include "Index.h"
// RAPsearch index
std::unordered_map<char, char> reducedIndexTable = {
    {'A', 'A'},
    {'C', 'C'},
    {'D', 'E'},
    {'E', 'E'},
    {'F', 'F'},
    {'G', 'G'},
    {'H', 'H'},
    {'I', 'I'},
    {'K', 'K'},
    {'L', 'I'},
    {'M', 'I'},
    {'N', 'E'},
    {'P', 'P'},
    {'Q', 'E'},
    {'R', 'K'},
    {'S', 'S'},
    {'T', 'S'},
    {'V', 'I'},
    {'W', 'F'},
    {'X', 'X'},
    {'Y', 'F'},
    {'Z', 'X'}};

Index::Index()
{
}

template <typename Out>
void split(const std::string &s, char delim, Out result)
{
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim))
    {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

std::string Index::AA2Reduced(std::string aa)
{
    std::string aaR;
    for (const char &i : aa)
    {
        aaR += reducedIndexTable[i];
    }
    return aaR;
}

void Index::training(std::string output)
{
    fasttext::FastText fasttext;
}

// typedef struct {
//     std::string key;
//     std::string value;
// } hash_index_table;

void Index::indexing(std::string finput, std::string output, int kmer, int label_index, bool isreduced)
{

    std::unordered_map<std::string, std::unordered_map<std::string, bool>> kmers;

    seqan::CharString id;
    seqan::CharString seq;
    std::ifstream input(finput);

    std::string label, prelabel, all_labels;
    std::string rProt; // protein with the reduced alphabet
    std::string kml;
    std::string fragment;

    int fragments = 0;

    std::ofstream fo(output + ".tr");
    // std::ofstream fos(output+".tri");
    std::string ks;
    int proteins = 0;
    seqan::SeqFileIn seqFileIn(seqan::toCString(finput));
    std::cout << "Processing input file ... " << std::endl;

    std::random_device rd;
    std::mt19937 rng(rd());
    int Sl;
    int rip;      // random position
    int l = 33;   // read length
    int k = kmer; // kmer size
    std::string read;

    while (!atEnd(seqFileIn))
    {
        std::cout << proteins << " processed reads"
                  << "\r";
        std::cout.flush();

        seqan::readRecord(id, seq, seqFileIn);

        if (isreduced)
        {
            rProt = AA2Reduced(seqan::toCString(seq));
        }
        else
        {
            rProt = seqan::toCString(seq);
        }

        prelabel = split(seqan::toCString(id), '|')[label_index];

        std::vector<std::string> _labels = split(prelabel, ';');
        std::stringstream _tagged_labels;

        for (int index_label = 0; index_label < _labels.size(); index_label++)
        {
            _tagged_labels << " __label__" << _labels[index_label] << "__ ";
        }

        label = _tagged_labels.str();

        Sl = rProt.length(); // length of the protein sequence
        proteins++;

        int min_kmers = int(Sl / k);
        int count_kmers = 0;
        std::uniform_int_distribution<int> uni(3, 5);

        for (int ix = 0; ix < k; ix++)
        {
            min_kmers = uni(rng);

            for (int i = ix; i <= Sl - k - 1; i += k)
            {
                ks = rProt.substr(i, k);
                kmers[ks][prelabel] = true;

                if (count_kmers % min_kmers == 0 && count_kmers > 0)
                {
                    fo << ks + ' ' << label << std::endl;
                    count_kmers = 0;
                }
                else
                {
                    fo << ks + ' ';
                    count_kmers++;
                }
            }

            count_kmers = 0;
        }
    }

    fo.close();

    std::ofstream fo2(output + ".kh");
    int ckmers = 0;
    for (const auto &arglabel : kmers)
    {

        fo2 << arglabel.first + "\t";
        for (const auto &lbl : arglabel.second)
        {
            fo2 << lbl.first << " ";
        }
        fo2 << std::endl;
        ckmers++;
    }

    std::cout << proteins << " proteins in the database " << std::endl;
    std::cout << ckmers << " unique " << kmer << "-mers " << std::endl;
    fo2.close();
}
