#include "tool.h"

double ticking()
{
  struct timeval timet;
  gettimeofday(&timet, NULL); 
  return timet.tv_sec +(timet.tv_usec/1000000.0);
}

uint64_t str_to_int (const char * str, size_t l)
{
  uint64_t strint = 0;
  for (size_t i = 0; i < l; i++) {
    uint8_t curr = 0;
    switch (str[i]) {
      case 'A': { curr = 0; break; }
      case 'T': { curr = 3; break; }
      case 'C': { curr = 1; break; }
      case 'G': { curr = 2; break; }
    }
    strint = strint << 2;
    strint = strint | curr;
  }
  return strint;
}

std::string int_to_str (uint64_t kmer, size_t l)
{
  uint8_t base;
  std::string str;
  for (int i=l; i>0; i--) {
    base = (kmer >> (i*2-2)) & 3ULL;
    char chr;
    switch(base) {
      case DNA_MAP::A: { chr = 'A'; break; }
      case DNA_MAP::T: { chr = 'T'; break; }
      case DNA_MAP::C: { chr = 'C'; break; }
      case DNA_MAP::G: { chr = 'G'; break; }
    }
    str.push_back(chr);
  }
  return str;
}


