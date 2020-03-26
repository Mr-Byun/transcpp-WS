/*********************************************************************************
*                                                                                *
*     fasta.cpp                                                                  *
*                                                                                *
*     contains methods for parsing a fasta file into a map of sequences          *
*                                                                                *
*********************************************************************************/

#include "fasta.h"

#include <boost/foreach.hpp>
//#include <boost/range/adaptor/map.hpp>
#include <limits>
#include <map>
#include <sstream>

#define foreach_ BOOST_FOREACH

typedef map<string, string>::iterator seqs_iter_type;

/*    Constructors    */

Fasta::Fasta() {}

Fasta::Fasta(string& fname) 
{
  file_name = fname;
  read(fname); 
}


/*    Getters   */

string& Fasta::getSeq(string& name)
{
  
  if (seqs.find(name) == seqs.end())
  {
    stringstream err;
    err << "ERROR: could not find seq with name " << name << endl;
    error(err.str());
    return seqs.end()->second;  // you will never get here!
  }
  else
    return seqs[name];
  
}

/*    I/O   */

void Fasta::read(string& fname)
{
  string line;
  ifstream file(fname.c_str());
  if (!file.good()) {
    stringstream err;
    err << "ERROR: Could not open file " << fname << endl;
    error(err.str());
  }
	
  string seq_name;
  
  if (file.is_open())
  {
    while(getline(file, line))
    {
      if (!line.empty())
      {
        if (line.at(0) == '>') // 
        {
          nseqs++;
          seq_name = line.substr(1,line.length()-1);
          names.push_back(seq_name);
          //cerr << seq_name << endl;
        } 
        else
        {
          seqs[seq_name] = seqs[seq_name] + line;
        }
      }
    }
  }
}

void Fasta::print(ostream& os)
{
  for(seqs_iter_type it = seqs.begin(); it != seqs.end(); it++)
  {
    
    os << '>' << it->first << endl;
    os << it->second << endl;
  }
}
          
          
          
