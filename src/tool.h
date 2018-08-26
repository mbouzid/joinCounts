#ifndef tools_h
#define tools_h

#include <sys/time.h>
#include <cstdlib>
#include <string>

enum DNA_MAP {A, C, G, T};  // A=1, C=0, T=2, G=3
static const char NUCLEOTIDES[4] = { 'A', 'C', 'G', 'T' };

double ticking();

uint64_t str_to_int (const char* str, size_t l);


std::string int_to_str(uint64_t kmer, size_t l);
#endif
