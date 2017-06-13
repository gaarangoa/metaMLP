#include "Signatures.h"
#include "functions.h"

using spp::sparse_hash_map;
using namespace std;

std::mutex mtx; 

// TODO!!: this split function is repeated in the INDEX.cpp, I have to move it to some common library so i can used it multiple times.
template<typename Out>
void splitx(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> splitx(const std::string &s, char delim) {
    std::vector<std::string> elems;
    splitx(s, delim, std::back_inserter(elems));
    return elems;
}
// END:TODO!!

// Constructor 
// Load the signatures from the input file
Signatures::Signatures(std::shared_ptr<fasttext::Args> a){
    
    args = a;

    file_model = a->smodel;
    fsignatures = file_model+".kh";
    kmer_size = a->kmer;
    seed_size = a->seed;
    isreduced = a->reduced;

    /////////////////////////////////////////////////////
    /*.............LOADING FASTTEXT MODEL.........*/
    /////////////////////////////////////////////////////
    
    std::cout<<"Loading Fast Text Model ... \n";
    fasttext.loadModel(file_model+".bin");
    
    //////////////////////////////////////////////////////
    /*.............LOADING SIGNATURES TO MEMORY.........*/
    //////////////////////////////////////////////////////
    std::cout<<"Loading Signatures ... \n";
    
    std::ifstream ifs(fsignatures);
    std::string line; 
    std::vector<std::string> iline;

    while( std::getline( ifs, line ) ){
        iline = splitx(line, '\t');
        master_signature_hash[iline[0].substr(0,seed_size)] = true;
        master_signature_hash_full[iline[0].substr(0,kmer_size)] = iline[1];
    }
    ifs.close();

}


void Signatures::predict(seqan::StringSet<seqan::Dna5String> &seqs, seqan::StringSet<seqan::CharString> &ids, std::vector<std::string>& readLabels, std::string& buffer, std::unordered_map < std::string, std::tuple < std::string, float > >& FuncPred){
    // TODO: This section has to be updated in the multiprocessing section'
    // Reads up to 10 records

    // fasttext::Vector queryVec(skipgram.args_->dim);
    // fasttext::Matrix wordVectors(skipgram.dict_->nwords(), skipgram.args_->dim);
    // skipgram.precomputeWordVectors(wordVectors);

    std::vector<std::string> readSeqs;
    std::vector<std::string> protSeqs;

    tsl::hopscotch_map< std::string, bool > signature_hash = master_signature_hash;
    tsl::hopscotch_map< std::string, std::string > signature_hash_full = master_signature_hash_full;

    int num_reads = 1;
    int total_reads = 0;
    
    int read_length = 100;

    /////////////////////////////////////////////////////
    /*.............GET READING FRAMES.........*/
    /////////////////////////////////////////////////////
    seqan::StringSet<seqan::String<seqan::AminoAcid>, seqan::Owner<seqan::ConcatDirect<> > > aaSeqs;
    if(isreduced){
        seqan::GeneticCode<seqan::MURPHY10> GCode; // reduce the aminoacid alphabet to 10
        seqan::translate(aaSeqs, seqs, seqan::SIX_FRAME, GCode);
    }else{
        seqan::GeneticCode<seqan::CANONICAL> GCode;
        seqan::translate(aaSeqs, seqs, seqan::SIX_FRAME, GCode);
    }
    

    // std::cout << aaSeqs[5] << std::endl;

    typedef seqan::Iterator< seqan::StringSet<seqan::String<seqan::AminoAcid>, seqan::Owner<seqan::ConcatDirect<> > > >::Type AIter;
    
    //////////////////////////////////////////////////////////////////////////////
    /*.............FILTER READS BY SIGNATURES AND CLASSIFY .........*/
    //////////////////////////////////////////////////////////////////////////////
    std::random_device rd;   
    std::mt19937 rng(rd());
    // std::uniform_int_distribution<int> uni(0, int(read_length/3)-seed_size-1);
    int frame = 1, rx;
    typedef seqan::Infix< seqan::String<seqan::AminoAcid> >::Type kmer;

    seqan::String<char> kimer;
    std::string toCSkmer;
    std::string KMER;
    std::string pre_buffer;

    std::size_t stop_c;

    std::string pkmer; 

    int l=0;
    int ishash;
    // mtx.lock();
    for( AIter it=begin(aaSeqs); it!=end(aaSeqs); ++it){

        KMER.clear();

        if(frame==6){
            frame=1;
            total_reads++;
        }else{
            frame++;
        }
        
        std::uniform_int_distribution<int> uni(0, length(*it)-seed_size-1);
        rx = uni(rng); // random position
        // Move iterator to a string-like structure
        
        seqan::move(kimer, *it);
        toCSkmer = seqan::toCString(kimer);
        
        // If the read has a stop codon, go to next reading frame:
        // stop_c = toCSkmer.find_first_of('*');
        if(stop_c<33) continue;
        l = toCSkmer.length();

        // Get kmer from a random position
        // KMER = toCSkmer.substr(rx, kmer_size);

        // ishash = signature_hash_full.count(KMER.substr(0,kmer_size));
        
        // make n tries to get the right kmer from the read
        int tries=0;
        // if(ishash>0){
            // ishash = 0;
            // rx=0;
            while(1){
                rx = uni(rng);
                KMER = toCSkmer.substr(rx, args->seed);
                ishash = signature_hash.count(KMER);
                if(ishash>0) break;
                if(tries == 20) break;
                tries++;
                // rx++;
            }
        // }

        // std::cout << ishash << "\t" << KMER << std::endl;

        // Got a kmer at all?, great, make a sentence and predict!! :)
        if(ishash>0){
            pre_buffer = KMER;
            int manykmers = 0;
            for(int ki=0; ki<20; ki++){
                rx = uni(rng);
                // TODO: Fixed kmer length for the substraction of subsequences. This parameter is fixed and is the same used for the training. 
                pkmer = toCSkmer.substr(rx, args->kmer);
                pre_buffer+=' '+pkmer;
                if(signature_hash_full.count(pkmer)>0){
                    manykmers++;
                }
                pkmer.clear();
            }

            if(manykmers>=args->mink){

                // buffer+=signature_hash_full[KMER]+' '+pre_buffer+'\n';
                buffer+=pre_buffer+'\n';
                pre_buffer.clear();
                
                readLabels.push_back(seqan::toCString(ids[total_reads]));
                std::stringstream iseq;
                iseq << seqs[total_reads];
                readSeqs.push_back(iseq.str());
                protSeqs.push_back(toCSkmer);
                num_reads++;
            }

        }

    }
    
    //  mtx.unlock();
    // std::unordered_map < std::string, std::tuple < std::string, float > > FuncPredLocal;
    std::stringstream trex(buffer);
    fasttext.predict(trex, 1, false, readLabels, 0, FuncPred, readSeqs);
    
    // std::cout << seqan::length(readLabels) << "\t" << seqan::length(readSeqs) << std::endl;
    // mtx.lock();
        // FuncPred = FuncPredLocal;
    // mtx.unlock();

    trex.str(std::string());
    readLabels.clear();
    buffer.clear();
    signature_hash.clear();
    seqan::clear(aaSeqs);
    
}



void Signatures::Display(std::string message){
    std::cout<<message<<"\n";
}

