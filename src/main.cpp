#include <cstdlib>
#include <cstdio> 
#include <zlib.h>

#include "kstring.h"
#include "kseq.h"
#include "tool.h"

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <exception>

#include <stdexcept>
#define kmer_length 31
#define buffer_size 100000

KSEQ_INIT(gzFile, gzread)

bool * term ;

struct s_stream
{ 
  char * fname;
  size_t id;
  gzFile  fp;
  kstream_t *ks;
};
typedef struct s_stream stream;




typedef std::map <uint64_t, uint32_t * > array;

struct s_matrix
{
  size_t n;
  array values;
};

typedef struct s_matrix matrix;


struct s_item
{
  uint64_t kmer;
  uint32_t count;
};

typedef struct s_item item;


int init (stream & s)
{

  s.fp = gzopen(s.fname, "r");
  if (not s.fp)
  { 
    std::cerr << "Failed to open "<< s.fname << std::endl; 
    exit(EXIT_FAILURE);
  }

  s.ks = ks_init(s.fp);

 
  return (EXIT_SUCCESS);
}


int read (stream & s, item & it) 
{
  int dret (0);
  kstring_t * str = (kstring_t*)calloc(1, sizeof(kstring_t));
  
  ks_getuntil(s.ks, 0, str, &dret);
  
  if ( not ks_eof(s.ks) )
  {
    
    std::string kmer (ks_release(str));
    uint64_t i_kmer = str_to_int(kmer.c_str(),kmer_length);
    
    ks_getuntil(s.ks, 0, str, &dret);
    uint32_t count = atoi(str->s);  

    it.kmer = i_kmer;
    it.count = count ;
    
    if (dret != '\n') while ((dret = ks_getc(s.ks)) > 0 && dret != '\n');
    dret = 1;
  }
  
  free(str->s);
  free(str);
  
  return dret;
  
}


void print (std::ostream & os, const matrix & counts, size_t a, size_t r)
{
  size_t n = counts.n;
  for (array::const_iterator i = counts.values.cbegin(); i != counts.values.cend(); ++i)
  {
    
    size_t nb (0);
    for (size_t j = 0; j < n; ++j)
    {
      if ((*i).second[j] >= a)
	++nb;

      if (nb >= r)
	break;
    }

    if (nb >= r)
    {
      os << int_to_str((*i).first,kmer_length) ;
      for (size_t j = 0; j < n; ++j)
	os << '\t' << (*i).second[j];
      os << '\n';
    }
  }
}

int not_terminated (const stream * s, size_t n)
{
  for (size_t i = 0; i < n; ++i)
    if (not term[i])
      return i;
  return -1;
}

int main (int argc, char * argv [])
{
  size_t r (1);
  size_t a (1); 
  long unsigned int b (1000);
  
  int c;
  while ((c = getopt(argc, argv, "r:a:b:")) >= 0) 
  {
    switch (c) 
    {
      case 'r': r = atoi(optarg); break;
      case 'a': a = atoi(optarg); break;
      case 'b': b = atoi(optarg); break;
    }
  }

  if (optind == argc) 
  {
    std::cerr << "Usage:\t" << argv[0] << " [options] <counts.tsv> [<counts.tsv> ...] " << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "\t-r INT\tmin recurrence [default:" << r << "]" << std::endl;
    std::cerr << "\t-a INT\tmin recurrence abundance [default:" << a << "]" << std::endl;
    std::cerr << "\t-b INT\tbuffer size [default:" << b << "]" << std::endl;

    exit (EXIT_FAILURE);
  }
  
  size_t n(argc-optind);

  stream * s = new stream [n];
  term = new bool[n];
  

  for (size_t i = 0; i < n; ++i)
  {
    s[i] = {argv[optind++],i,NULL,NULL};
    init(s[i]);
    term[i] = false; 
  }
  
  matrix counts = {n, std::map <uint64_t,uint32_t *>()}; 
  matrix last = {n,std::map<uint64_t, uint32_t *>()};
  int m (0);

  while (( m = not_terminated(s,n) )!= -1)
  {
    size_t j (0);
    item max;
  
    while ( (not term[m]) and (j < b))
    {
      int stat = read(s[m],max);
      
      if (not stat)
      {
	term[m] = true;
	break;
      }

      if (counts.values.find(max.kmer) == counts.values.end())
      {
	counts.values.emplace(max.kmer,new uint32_t[counts.n]);
      }

      counts.values[max.kmer][m] = max.count;

      ++j;
    }

    std::map<uint64_t,uint32_t * >::iterator p = last.values.begin();
    if ( not last.values.empty() )
    {
      while ( (p != last.values.end()) and ((*p).first < max.kmer))
      {
	if (counts.values.find((*p).first) == counts.values.end())
	{
	  counts.values.emplace((*p).first,new uint32_t[counts.n]);
	}
    
	for (size_t o = 0; o < n; ++o)
	{
	  if ( (*p).second[o] != 0)
	    counts.values[(*p).first][o] = (*p).second[o];
	}

	p = last.values.erase(p);
      }
    }

    
       
    for (int i = 0; i < (int)n; ++i)
    {
      
      if ((not term[i]) and (i!=m))
      {
	item next;	
	
	int status ;
		
	while ( not term[i] ) 
	{
	  status = read(s[i],next);
 
	  if (not status)
	  {
	    term[i]= true;
	    break;
	  }

	    


 
	  if (next.kmer > max.kmer)
	  {
	       
	    break;
	  }
	  
	 

	  
	  if (counts.values.find(next.kmer) == counts.values.end())
	  {
	    counts.values.emplace(next.kmer,new uint32_t[counts.n]);
	  }
	  
	  counts.values[next.kmer][i] = next.count; 
	   
	}

	if (status)
	{
	  if (last.values.find(next.kmer) == last.values.end())
	  {
	    last.values.emplace(next.kmer,new uint32_t[counts.n]);
	  }
	  
	  last.values[next.kmer][i] = next.count; 
	   
	   
	}

      }

    }

    print(std::cout,counts,a,r);
    counts.values.clear();
     
  }

  for (size_t i = 0; i < n; ++i)
    ks_destroy (s[i].ks);
  
  delete [] s;
  delete [] term;
  return (EXIT_SUCCESS);
}
