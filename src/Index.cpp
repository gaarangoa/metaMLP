#include "Index.h"
// RAPsearch index
std::unordered_map<char, char> reducedIndexTable = {
                                                        {'A','A'},
                                                        {'C','C'},
                                                        {'D','E'},
                                                        {'E','E'},
                                                        {'F','F'},
                                                        {'G','G'},
                                                        {'H','H'},
                                                        {'I','I'},
                                                        {'K','K'},
                                                        {'L','I'},
                                                        {'M','I'},
                                                        {'N','E'},
                                                        {'P','P'},
                                                        {'Q','E'},
                                                        {'R','K'},
                                                        {'S','S'},
                                                        {'T','S'},
                                                        {'V','I'},
                                                        {'W','F'},
                                                        {'X','X'},
                                                        {'Y','F'},
                                                        {'Z','X'}
                                                    };

// DIAMOND index
// std::unordered_map<char, char> reducedIndexTable = {
//                                                         {'A','A'},
//                                                         {'C','C'},
//                                                         {'D','K'},
//                                                         {'E','K'},
//                                                         {'F','F'},
//                                                         {'G','G'},
//                                                         {'H','H'},
//                                                         {'I','I'},
//                                                         {'K','K'},
//                                                         {'L','I'},
//                                                         {'M','M'},
//                                                         {'N','K'},
//                                                         {'P','P'},
//                                                         {'Q','K'},
//                                                         {'R','K'},
//                                                         {'S','A'},
//                                                         {'T','A'},
//                                                         {'V','I'},
//                                                         {'W','W'},
//                                                         {'X','X'},
//                                                         {'Y','Y'},
//                                                         {'Z','X'}
//                                                     };


Index::Index(){
    
}

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

std::string Index::AA2Reduced(std::string aa){
    std::string aaR;
    for(const char& i: aa){
        aaR+=reducedIndexTable[i];
    }
    return aaR;
}

void Index::training(std::string output){
    fasttext::FastText fasttext; 
}


// typedef struct {
//     std::string key;
//     std::string value;
// } hash_index_table;


void Index::indexing(std::string finput, std::string output, int kmer, int label_index, bool isreduced){
    
    std::unordered_map< std::string,std::unordered_map< std::string, bool > > kmers;

    seqan::CharString id;
    seqan::CharString seq;
    std::ifstream input(finput);
    
    std::string label, prelabel;
    std::string rProt; // protein with the reduced alphabet
    std::string kml;
    std::string fragment;

    int fragments=0;

    std::ofstream fo(output+".tr");
    // std::ofstream fos(output+".tri");
    std::string ks;
    int proteins=0;
    seqan::SeqFileIn seqFileIn(seqan::toCString(finput));
    std::cout << "Processing input file ... " << std::endl;


    std::random_device rd;   
    std::mt19937 rng(rd());
    int Sl;
    int rip; // random position
    int l = 33; // read length
    int k = kmer; // kmer size
    std::string read;
    
    while (!atEnd(seqFileIn)){
        std::cout << proteins << " processed reads" <<"\r";
        std::cout.flush();

        seqan::readRecord(id, seq, seqFileIn);
        
        if(isreduced){
            rProt =  AA2Reduced(seqan::toCString(seq));
        }else{
            rProt = seqan::toCString(seq);}
        
        prelabel = split(seqan::toCString(id), '|')[label_index];
        label = "__label__"+prelabel+"__";
        
        // TODO: the i+=2 takes each protein and slides the window with two amnoacids. This parameter is set to 2 to avoid to get too many "reads" that are used for training. 

        Sl = rProt.length(); // length of the protein sequence
        // for(int i=0; i<Sl-l; i+=3){

        //     for (int ri=0; ri<10; ri++){
        //         std::uniform_int_distribution<int> uni(i, i+l-k);
        //         rip = uni(rng);
        //         ks = rProt.substr(rip, k);
        //         kmers[ks][prelabel]=true;
        //         read += ' '+ks;
        //     }

        //     // fo << label << "\t" << prelabel+' ' << read << std::endl;
        //     // fo << label << "\t" << read << std::endl;
        //     read.clear();

        // }

        // read.clear();
        // fragments = 0;
        proteins++;

        // for computing the skipgram I follow a different approach, first  I don't need to split the sequence into reads, and just take the kmers and store them into a file that will be used for the training of the skipgram. 
        
        // TODO replace this if you want to put in index the label information
        // fo << label << "\t" << prelabel+' ';

        fo << label << '\t';

        for (int i=0; i<=Sl-k; i++){
            ks = rProt.substr(i, k);
            kmers[ks][prelabel]=true;
            // fos << ks+' ';
            if(i % 20 == 0 && i>0){
                fo << ks+' ' << std::endl;
                fo << label << '\t';
                // TODO: here add it if you want to put the labels information
                // fo << label << "\t" << prelabel+' ';
            }else{
                fo << ks+' ';
            }
        }

        fo << std::endl;
        // fos << std::endl;
        
    }

    fo.close();
    // fos.close();

    // Store kmers to a text file
    // hash_index_table HASH_INDEX[kmers.size()];

    std::ofstream fo2(output+".kh");
    int ckmers=0;
    for(const auto& arglabel: kmers){

        fo2 << arglabel.first + "\t";
        for(const auto& lbl: arglabel.second){
            fo2 << lbl.first << " " ;
        }
        fo2 << std::endl;

        // HASH_INDEX[ckmers].key = arglabel.first;
        // HASH_INDEX[ckmers].value = arglabel.first;

        ckmers++;
    }

    std::cout << proteins << " proteins in the database " << std::endl;
    std::cout << ckmers << " unique " << kmer << "-mers " << std::endl;
    fo2.close();
    
    // int length_hash_index = sizeof(HASH_INDEX)/sizeof(hash_index_table);
    // hsize_t dim[1];
    // dim[0] = sizeof(HASH_INDEX)/sizeof(hash_index_table);
    // int rank = sizeof(dim)/sizeof(hsize_t);

    // // defining data to pass to HDF5
    // H5::CompType mtype(sizeof(hash_index_table));
    
    // mtype.insertMember("key", HOFFSET(hash_index_table, key), H5::StrType());
    // mtype.insertMember("value", HOFFSET(hash_index_table, value), H5::StrType());

    // // dataspace preparation
    // H5::DataSpace space(rank, dim);
    // H5::H5File *file = new H5::H5File(output+".kh.hdf", H5F_ACC_TRUNC);
    // H5::DataSet *dataset = new H5::DataSet(file->createDataSet("kmer_index", mtype, space));

    // delete dataset;
    // delete file;



}
