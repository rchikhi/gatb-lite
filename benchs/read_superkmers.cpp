
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

typedef vector<shared_ptr<vector<string>>> read_packets_t ;

void chop_read_into_kmers(string& read)
{
    int k = 100;
    if (read.size() < k) return;

    for (int i = 0; i < read.size()-k+1; i++)
    {
        std::cout << read.substr(i,i+k) << std::endl;
    }
}

void reads_to_superkmer(read_packets_t &read_packets, std::mutex& lock, int thread_id)
{
    // critical section
    lock.lock();
    auto read_packet = std::move(read_packets.back());
    read_packets.pop_back();
    lock.unlock();

    for (auto read: *read_packet)
    {
        //std::cout << read << std::endl;
        chop_read_into_kmers(read);
    }
    
    read_packet.reset(); // would like to do it but somehow cant as it segfaults
}

int
main(int argc, char** argv)
{
    if (argc < 2)
        return 1;

    int nb_threads=6;
    ThreadPool pool(nb_threads);
    shared_ptr<vector<string>> read_packet = make_shared<vector<string>>();
    read_packets_t read_packets;
    std::mutex lock;

    process_fastq(argv[1], [&read_packets, &read_packet, &pool, &lock](fastq_record<>& rec) { 

            string&& read_str = string(rec.sequence().begin(),rec.sequence().end());
            read_packet->push_back(std::move(read_str));

            if (read_packet->size() == 100)
            {
               auto r_t_s_wrapper = [&read_packets, &lock] (int thread_id) 
               { reads_to_superkmer(read_packets, lock, thread_id);};

               lock.lock();
               read_packets.push_back(std::move(read_packet));
               lock.unlock();

               read_packet = make_shared<vector<string>>();
               
               pool.enqueue(r_t_s_wrapper);
               //r_t_s_wrapper(0); // single threaded
            }
		  });
    pool.join();

}
