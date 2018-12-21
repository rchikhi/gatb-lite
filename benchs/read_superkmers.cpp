
#include <cstddef>
#include <iostream>
#include <cstring>

#include "gatbl/sys/file.hpp"
#include "gatbl/fastx.hpp"

// superkmer-specific stuff
#include "ThreadPool.h"

using namespace gatbl;

template<typename F>
inline void
process_fastq(std::string filename, F&& callback)
{
    sys::file_descriptor fd(filename);
    auto data = fd.mmap<const char>();
    // FIXME: those are optional but we really want to knwow if they work for testing
    assert(data.advise_hugepage());
    assert(data.advise_sequential());
    for (auto& rec : seq_record_subrange<fastq_record<const char*>>(data)) {
        callback(rec);
    }
}

using namespace std;

typedef shared_ptr<vector<string>> read_packet_t ;
typedef vector<read_packet_t> read_packets_t ;
typedef uint64_t kmer_type;

inline static char char_to_nt(char c)
{
    return (c>>1) & 3;
}

inline static char char_to_nt_rc(char c)
{
    return ((c>>1)^2) & 3;
}

// to avoid being optimized out
template<typename T>
void
use(T&& t)
{
        __asm__ __volatile__("" ::"g"(t));
}


void emit_kmer(kmer_type kmer)
{
    // do what we want with that kmer
    //std::cout << kmer << std::endl;
    use(kmer);
}

bool is_minimizer_allowed(uint32_t mmer)
{
    // TODO
    // for now is strict lexical order, not KMC's order
    return true;
}

inline static uint64_t revcomp (const uint64_t & x, int k)
{
    uint64_t res = x;
    res = ((res>> 2 & 0x3333333333333333) | (res & 0x3333333333333333) <<  2);
    res = ((res>> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) <<  4);
    res = ((res>> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) <<  8);
    res = ((res>>16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res>>32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
    res = res ^ 0xAAAAAAAAAAAAAAAA;

    return (res >> (2*(32-k))) ;
}

void init_minimizers(int _minimizerSize, vector<uint32_t> &_mmer_lut)
{
    uint32_t _mask  = (1 << (2*_minimizerSize)) - 1;

    uint32_t nbminims_total = (1 << (2*_minimizerSize));

    for(uint32_t ii=0; ii< nbminims_total; ii++)
    {
        uint32_t mmer = ii;
        uint32_t rev_mmer = revcomp(mmer, _minimizerSize);
        if(rev_mmer < mmer) mmer = rev_mmer;

        if (!is_minimizer_allowed(mmer),_minimizerSize)
            mmer = _mask;

        _mmer_lut[ii] = mmer;
    }
}

void compute_minimizer(uint64_t kmer, uint32_t &minimizer, unsigned int& minimizer_position, unsigned int k, vector<uint32_t> &_mmer_lut, unsigned int _minimizerSize)
{
    minimizer = 0;
    uint32_t mmer;
    uint32_t  mmerMask  = (1 << (2*_minimizerSize)) - 1;
    for (size_t i=0; i<k-_minimizerSize+1; ++i)
    {
        mmer    = _mmer_lut[(kmer >>i)& mmerMask];
        if (mmer < minimizer)
        {
            minimizer = mmer;
            minimizer_position = i;
        }
    }
}

void chop_read_into_kmers(string& read, vector<uint32_t> &_mmer_lut, unsigned int k, unsigned int _minimizerSize, vector<int> &superkmer_positions, vector<unsigned int> &nb_kmers_thread, int thread_id)
{
    kmer_type kmer = 0, kmer_rc = 0;
    kmer_type kmerMask = (1LL << (k*2)) - 1;
    uint32_t  mmerMask  = (1 << (2*_minimizerSize)) - 1;
    unsigned int minimizer_position;
    uint32_t minimizer, mmer;

    // TODO skip kmers contaning N's, for now theyre transformed into G's
    
    for (size_t i=0; i<k; ++i)
    {
        char c  = char_to_nt(read[i]);
        kmer    = (kmer<<2) + c;
        kmer_rc = (kmer>>2) + (c^2);
    }
    emit_kmer(std::min(kmer,kmer_rc));
    nb_kmers_thread[thread_id]++;
    compute_minimizer(kmer, minimizer, minimizer_position, k, _mmer_lut, _minimizerSize);

    size_t read_size = read.size();
    for (size_t i = k; i < read_size; i++)
    {
        char c  = char_to_nt(read[i]);
        kmer    = (( kmer << 2) +  c) & kmerMask;
        kmer_rc = (( kmer >> 2) +  (c^2)) & kmerMask;
        emit_kmer(std::min(kmer,kmer_rc));
        nb_kmers_thread[thread_id]++;

        mmer    = _mmer_lut[kmer & mmerMask];
        minimizer_position--;
        if (mmer <= minimizer)
        {
            if (mmer != minimizer) 
            {
                superkmer_positions.push_back(i);
                minimizer = mmer;
            }
            minimizer_position = i-_minimizerSize+1;
        }
        else
        {
            if (minimizer_position < 0)
            {
                uint32_t old_minimizer = minimizer;
                compute_minimizer(kmer, minimizer, minimizer_position, k, _mmer_lut, _minimizerSize);
                if (minimizer != old_minimizer)
                    superkmer_positions.push_back(i);
            }
        }
    }
}

void emit_superkmers(string &read, vector<int> &superkmer_positions, int k, vector<unsigned int> &nb_superkmers_thread, int thread_id)
{
    // do what we want with that kmer
    //std::cout << kmer << std::endl;
    
    size_t read_size = read.size();
    string superkmer = "";
    unsigned int next_break_i = 0;
    unsigned int next_break;
    if (superkmer_positions.size() == 0)
       next_break = read.size()+1;
    else
       next_break = superkmer_positions[next_break_i];
    uint64_t last_kmer = 0;
    kmer_type kmerMask = (1LL << (k*2)) - 1;

    for (size_t i = 0; i < read_size; i++)
    {
        char c  = char_to_nt(read[i]);
        last_kmer    = ((last_kmer<<2) + c ) & kmerMask;
        superkmer += c;

        if (i == next_break)
        {
            use(superkmer); // emit the superkmer here
            nb_superkmers_thread[thread_id]++;

            next_break_i++;
            next_break = superkmer_positions[next_break_i];
            superkmer = last_kmer; // it works?! TODO check
        }
    }
    
    use(superkmer); // emit the last superkmer
    nb_superkmers_thread[thread_id]++;
}


void reads_to_superkmer(read_packet_t read_packet, 
        vector<uint32_t> _mmer_lut, unsigned int k, unsigned int _minimizerSize, int thread_id,
        vector<unsigned int> &nb_kmers_thread, vector<unsigned int> &nb_superkmers_thread)
{
    for (auto read: *read_packet)
    {
        if (read.size() < k) continue;

        //std::cout << read << std::endl;
        vector<int> superkmer_positions;
        chop_read_into_kmers(read, _mmer_lut, k, _minimizerSize, superkmer_positions, nb_kmers_thread, thread_id);
        emit_superkmers(read, superkmer_positions, k, nb_superkmers_thread, thread_id);
    }
    
    read_packet.reset(); 
}

int
main(int argc, char** argv)
{
    if (argc < 2)
        return 1;

    int nb_threads = 6;

    if (argc == 3)
        nb_threads = atoi(argv[2]);

    int _minimizerSize = 10;
    int k = 25;

    ThreadPool pool(nb_threads);
    shared_ptr<vector<string>> read_packet = make_shared<vector<string>>();
    read_packet->reserve(100000);
    vector<unsigned int> nb_kmers_thread(nb_threads), nb_superkmers_thread(nb_threads);
    uint64_t nb_kmers = 0, nb_superkmers = 0;
    uint32_t nbminims_total = (1 << (2*_minimizerSize));
    vector<uint32_t> _mmer_lut(nbminims_total); 
    init_minimizers(_minimizerSize, _mmer_lut);

    for (int i = 0; i < nb_threads; i++)
    { nb_kmers_thread[i] = 0; nb_superkmers_thread[i] = 0;}

    process_fastq(argv[1], [&read_packet, &pool, &nb_kmers_thread, &nb_superkmers_thread, &_mmer_lut, k, _minimizerSize](fastq_record<>& rec) { 

            string&& read_str = string(rec.sequence().begin(),rec.sequence().end());
            read_packet->push_back(std::move(read_str));

            if (read_packet->size() == 100000)
            {
               auto r_t_s_wrapper = [read_packet, &nb_kmers_thread, &nb_superkmers_thread, _mmer_lut, k, _minimizerSize] (int thread_id) 
               { reads_to_superkmer(read_packet, _mmer_lut, k, _minimizerSize, thread_id, nb_kmers_thread, nb_superkmers_thread);};

               
               pool.enqueue(r_t_s_wrapper);
               //r_t_s_wrapper(0); // single threaded
               
               read_packet = make_shared<vector<string>>();
               read_packet->reserve(100000);
            }
		  });
    pool.join();

    for (int i = 0; i < nb_threads; i++)
    { nb_kmers += nb_kmers_thread[i] ; nb_superkmers += nb_superkmers_thread[i];}

    std::cout << "nb kmers:      " << nb_kmers      << std::endl;
    std::cout << "nb superkmers: " << nb_superkmers << std::endl;
}
