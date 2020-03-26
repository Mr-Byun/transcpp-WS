/*********************************************************************************
*                                                                                *
*     twobit.h                                                                   *
*                                                                                *
*     Contains methods to read sequences from .2bit compressed genomes           *
*                                                                                *
*********************************************************************************/

# ifndef TWOBIT_H
# define TWOBIT_H

# include <iostream>
# include <fstream>
# include <cstdlib>
# include <string>
# include <vector>
# include <boost/cstdint.hpp>
# include <cmath>
# include <map>

using namespace std;

struct Header
{
  unsigned long signature;
  unsigned long version;
  unsigned long sequenceCount;
  unsigned long reserved;
};


struct IdxEntry
{
  int nameSize;
  string name;
  unsigned long offset;
};


struct Record
{
  string name;
  int dnaSize;
  int nBlockCount;
  vector<int> nBlockStarts;
  vector<int> nBlockSizes;
  int maskBlockCount;
  vector<int> maskBlockStarts;
  vector<int> maskBlockSizes;
  int reserved;
  unsigned long offset;
};
  
class TwoBit
{
private:
  string filename;
  Header header;
  vector<IdxEntry> index; 
  vector<Record>   records;
  
  void readHeader(void);
  void readIndex(void);
  void readRecords(void);
  
  void byte2seq(char, char*);
  int getN(string & chr); 
  
public:
  // Constructors
  TwoBit();
  TwoBit(string& fname);
  
  void read(string& fname);
  
  // Print
  void printHeader(ostream & os);
  void printIndex(ostream & os);
  
  // Getters
  string getSequence(string & chr);
  string getSequence(string & chr, int start, int end);
  string getFileName() { return filename; }
  
  
};


#endif
