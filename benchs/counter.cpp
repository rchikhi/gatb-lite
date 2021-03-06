
#include <cstddef>
#include <iostream>
#include <cstring>

#include "gatbl/sys/file.hpp"
#include "gatbl/fastx.hpp"

using namespace gatbl;

template<typename F>
inline void
process_fastq(std::string filename, F&& callback)
{
    sys::file_descriptor fd(filename);
    auto data = fd.mmap<const char>();
    // FIXME: those are optional but we really want to knwow if they work for testing
    data.advise_hugepage();
    data.advise_sequential();
    for (auto& rec : seq_record_subrange<fastq_record<const char*>>(data)) {
        callback(rec);
    }
}

template<typename T>
void
use(T&& t)
{
    __asm__ __volatile__("" ::"g"(t));
}

void
clobber()
{
    asm volatile("" : : : "memory");
}

using namespace std;

int
main(int argc, char** argv)
{
    if (argc < 2)
        return 1;
    size_t count = 0;
    process_fastq(argv[1], [&](fastq_record<>& rec) { 
		    use(rec);
		    clobber();
            count++;
		  });

    std::cout << count << std::endl;
    return 0;
}
