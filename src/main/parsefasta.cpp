/*********************************************************************************
*                                                                                *
*     parse_fasta.cpp                                                            *
*                                                                                *
*     main file for parsing fast sequences                                       *
*                                                                                *
*********************************************************************************/

#include "fasta.h"

#include <unistd.h>
#include <getopt.h>

static const char *optString = "hi:p:";

static const struct option longOpts[] = {
    { "help",        no_argument,       NULL, 'h' },
    { "input-file",  required_argument, NULL, 'i' },
    { "print",       no_argument,       NULL, 'p' }
};

void display_usage()
{
  cerr << endl << "\t Usage" << endl << endl
       << "\t parsefasta [options] -i [infile]" << endl << endl
       << "\t Options" << endl
       << "\t --print    [-p]   print the output file" << endl << endl;
  exit(1);
}



int main(int argc, char* argv[])
{
  int opt = 0;
  int longIndex = 0;
  
  bool print = false;
  
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
      case 0: 
        
        break;
    }
    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  }
  
  Fasta in(infile_name);
  
  if (print)
   in.print(cout);
}
  
