/*********************************************************************************
*                                                                                *
*     getsequence.cpp                                                            *
*                                                                                *
*     gets a sequence from a 2bit file using coordinates                         *
*                                                                                *
*********************************************************************************/

#include "twobit.h"

#include <unistd.h>
#include <getopt.h>

static const char *optString = "hi:pc:s:e:";

static const struct option longOpts[] = {
    { "help",        no_argument,       NULL, 'h' },
    { "input-file",  required_argument, NULL, 'i' },
    { "print",       no_argument,       NULL, 'p' },
    { "chr",         required_argument, NULL, 'c' },
    { "start",       required_argument, NULL, 's' },
    { "end",         required_argument, NULL, 'e' }
};

void display_usage()
{
  cerr << endl << "\t Usage" << endl << endl
       << "\t getseq [options] -i [infile]" << endl << endl
       << "\t Options" << endl
       << "\t --print    [-p]   print the header for the file" << endl
       << "\t --chr      [-c]   the chromosome or header to read from" << endl
       << "\t --start    [-s]   the starting base" << endl
       << "\t --end      [-e]   the ending base" << endl << endl;
  exit(1);
}

int main(int argc, char* argv[])
{
  int opt = 0;
  int longIndex = 0;
  
  bool   print = false;
  string chr;
  int    start;
  int    end;
  
  string infile_name;
  
  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while(opt != -1)
  {
    switch (opt)
    {
      case 'h':
        display_usage();
        break;
      case 'i':
        infile_name = optarg;
        break;
      case 'p':
        print = true;
        break;
      case 'c':
        chr = optarg;
        break;
      case 's':
        start = atoi(optarg);
        break;
      case 'e':
        end = atoi(optarg);
        break;
      case 0: 
        
        break;
    }
    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  }
  
 
  TwoBit genome(infile_name);
  
  if (print)
  {
    genome.printHeader(cout);
    genome.printIndex(cout);
    return(0);
  }
 
  cout << genome.getSequence(chr, start, end) << endl;
 
  return 0;
}
