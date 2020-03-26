/*********************************************************************************
*                                                                                *
*     fasta.h                                                                    *
*                                                                                *
*     contains methods for parsing a fasta file into a map of sequences          *
*                                                                                *
*********************************************************************************/

#ifndef FASTA_H
#define FASTA_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <boost/cstdint.hpp>
#include <cmath>
#include <map>
#include <boost/unordered_map.hpp>

#include "utils.h"

using namespace std;

class Fasta 
{
private:
  string file_name;
  int    nseqs;
  vector<string>      names;
  map<string, string> seqs;
  
public:
  // Constructors
  Fasta();
  Fasta(string& fname);
  
  // Getters
  string& getSeq(string& name);
  string& getFileName() { return file_name; }
  vector<string>& getNames() { return names; }
  
  // I/O
  void read(string& fname);
  void write(string& fname);
  void print(ostream& os);

};

#endif
