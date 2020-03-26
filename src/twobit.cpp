/*********************************************************************************
*                                                                                *
*     twobit.cpp                                                                 *
*                                                                                *
*     Contains methods to read sequences from .2bit compressed genomes           *
*                                                                                *
*********************************************************************************/

#include "twobit.h"
#include "utils.h"
#include <sstream>

using namespace std;
    
/* simple utility function mapping integer to char */
char seqmap(int i)
{
  if (i == 0)
    return 'T';
  else if (i == 1)
    return 'C';
  else if (i == 2)
    return 'A';
  else
    return 'G';
}

/*  Constructor   */

TwoBit::TwoBit() {}

TwoBit::TwoBit(string& fname)
{
  read(fname);
}

void TwoBit::read(string& fname)
{
  filename = fname;
  readHeader();
  readIndex();
  readRecords();
}
  
/*    Read File format    */

void TwoBit::readHeader(void)
{
  ifstream myFile(filename.c_str(), ios::in | ios::binary);
  uint32_t x;
  
  stringstream err;
  if (!myFile.read((char*)&x, sizeof(uint32_t)))
  {
    err << "ERROR: readHeader() could not read file " << filename << endl;
    error(err.str());
  }
  if (x != 0x1A412743)
  {
    err << "ERROR: readHeader() architecture did not return 0x1A412743" << endl;
    error(err.str());
  }
  header.signature = x;
  myFile.read((char*)&x, sizeof(uint32_t));
  if (x != 0 ) 
  {
    err << "ERROR: readHeader() version (" << x << ") not equal to 0!" << endl;
    error(err.str());
  }
  header.version = x;
  myFile.read((char*)&x, sizeof(uint32_t));
  header.sequenceCount = x;
  myFile.read((char*)&x, sizeof(uint32_t));
  header.reserved = x;
  myFile.close();
}

void TwoBit::readIndex(void)
{
  ifstream myFile(filename.c_str(), ios::in | ios::binary);
  myFile.seekg(16);
  
  unsigned int seqCount = header.sequenceCount;
  index.resize(seqCount);
  for (unsigned int i=0; i<seqCount; i++)
  {
    uint8_t x;
    myFile.read((char*)&x, sizeof(uint8_t));
    index[i].nameSize = (int) x;
    
    char * buffer = new char[index[i].nameSize];
    myFile.read(buffer, index[i].nameSize);
    index[i].name.assign(buffer, index[i].nameSize);
    delete[] buffer;
    
    uint32_t off;
    myFile.read((char*)&off, sizeof(uint32_t));
    index[i].offset = off;
  }
  myFile.close();
}

void TwoBit::readRecords(void)
{
  ifstream myFile(filename.c_str(), ios::in | ios::binary);
  records.resize(index.size());
  uint32_t x;
  
  unsigned int size = index.size();
  for (unsigned int i=0; i<size; i++)
  {
    records[i].name = index[i].name;
    myFile.seekg(index[i].offset);
    
    myFile.read((char*)&x, sizeof(uint32_t));
    records[i].dnaSize = x;

    myFile.read((char*)&x, sizeof(uint32_t));
    records[i].nBlockCount = x;
    
    records[i].nBlockStarts.reserve(records[i].nBlockCount);
    records[i].nBlockSizes.reserve(records[i].nBlockCount);
    
    for (int j=0; j<records[i].nBlockCount; j++)
    {
      myFile.read((char*)&x, sizeof(uint32_t));
      records[i].nBlockStarts.push_back(x);
    }
    for (int j=0; j<records[i].nBlockCount; j++)
    {
      myFile.read((char*)&x, sizeof(uint32_t));
      records[i].nBlockSizes.push_back(x);
    }
    
    myFile.read((char*)&x, sizeof(uint32_t));
    records[i].maskBlockCount = x;
    
    records[i].maskBlockStarts.reserve(records[i].maskBlockCount);
    records[i].maskBlockSizes.reserve(records[i].maskBlockCount);
    
    for (int j=0; j<records[i].maskBlockCount; j++)
    {
      myFile.read((char*)&x, sizeof(uint32_t));
      records[i].maskBlockStarts.push_back(x);
    }
    for (int j=0; j<records[i].maskBlockCount; j++)
    {
      myFile.read((char*)&x, sizeof(uint32_t));
      records[i].maskBlockSizes.push_back(x);
    }
    
    myFile.read((char*)&x, sizeof(uint32_t));
    records[i].reserved = 0;
    
    records[i].offset = myFile.tellg();;
  }
}
    

int TwoBit::getN(string & chr)
{
  unsigned int size = index.size();
  for (unsigned int i=0; i<size; i++)
  {
    if (records[i].name == chr)
      return(i);
  }
  stringstream err;
  err << "ERROR: getN() could not find chromosome " << chr << " in genome" << endl;
  error(err.str());
  return 0; // you will not get here!
}


void TwoBit::byte2seq(char x, char* out)
{
  int c = (x >> 4) & 0xf;
  int f = x & 0xf;
  int cc = (c >> 2) & 0x3;
  int cf = c & 0x3;
  int fc = (f >> 2) & 0x3;
  int ff = f & 0x3;

  out[0]=seqmap(cc);
  out[1]=seqmap(cf);
  out[2]=seqmap(fc);
  out[3]=seqmap(ff);
}

// return entire chromosome
string TwoBit::getSequence(string& chr)
{
  int idx = getN(chr); // get index of chromosome
  int start = 0;
  int end   = records[idx].dnaSize;
  
  return getSequence(chr, start, end);

}
  
string TwoBit::getSequence(string& chr, int start, int end)
{
  int idx = getN(chr); // get index of chromosome
  
  // error checking
  if (start > end)
  {
    cerr << "WARNING: getSequence() start position is greater than end position." << endl
         << "         Reversing start and end coordinates." << endl;
    int tstart = start;
    start = end;
    end = tstart;
  }
  if ( (start < 0) | (end < 0) )
  {
    stringstream err;
    err << "ERROR: getSequence() attempted to return negative coordinates" << endl;
    error(err.str());
  }
  if ( (end > records[idx].dnaSize) | (start > records[idx].dnaSize) )
  {
    stringstream err;
    err << "ERROR: getSequence() attempted to return sequence longer than chromosome" << endl;
    error(err.str());
  }
  
  unsigned long length = end - start + 1;
  unsigned long fstart; // file starting read position 
  unsigned long flen;
  unsigned long fend;
  int trim;
        
  fstart = (start-1)/4 + records[idx].offset; // note that counting starts at 0, so we subtract 1
  fend = fstart + (unsigned long) ceil(length/4.0);
  flen = fend - fstart;
  trim = (start-1)%4; // the characters at the beginning to trim off
  ifstream myFile(filename.c_str(), ios::in | ios::binary);
  myFile.seekg(fstart);
  
  char * buffer = new char[flen];
  char * curbyte = new char[4];
  string out;
  out.resize(flen*4);
  
  myFile.read(buffer, flen);
  //cerr << "sequence is " << length << " bases long" << endl;
  //cerr << "started reading at char " << fstart << ". Will need to trim " << trim << " from beginning" << endl;
  for (unsigned long i = 0; i<flen; i++)
  {
    byte2seq(buffer[i], curbyte);
    //out.append(curbyte,4);
    //cerr << i*4 << "\t" << curbyte[0] << endl;
    out[i*4]     = curbyte[0];
    //cerr << i*4 + 1 << "\t" << curbyte[1] << endl;
    out[i*4 + 1] = curbyte[1];
    //cerr << i*4 + 2  << "\t" << curbyte[2] << endl;
    out[i*4 + 2] = curbyte[2];
    //cerr << i*4 + 3 << "\t" << curbyte[3] << endl;
    out[i*4 + 3] = curbyte[3];
    
  }
  delete[] buffer;
  delete[] curbyte;
  
  myFile.close();
  
  //cerr << out << endl;
  
  out = out.substr(trim, length);
  
  //cerr << "trimmed to " << out << endl;
  
  /* now I need to figure out if we are over a stretch of Ns and replace
  the out sequence accordingly */
  
  int sublen, last, substart, subend;
  for (int i = 0; i<records[idx].nBlockCount; i++)
  {
    // check to see if the block starts before this
    if (records[idx].nBlockStarts[i] < start)
    {
      last = records[idx].nBlockStarts[i]+records[idx].nBlockSizes[i];
      // check to see if this extends into the retrieved region
      if ( last >= start)
      {
        substart = 0;
        sublen = min(last, end)-start+1;
        //cerr << "replacing bases " << substart << " through " << sublen << " with N" << endl;
        out.replace(substart,sublen,sublen,'N');
        //cerr << out << endl;
      }
    } else if (records[idx].nBlockStarts[i] <= end) 
    {
      last = records[idx].nBlockStarts[i]+records[idx].nBlockSizes[i];
      substart = records[idx].nBlockStarts[i]-start+1;
      subend = min(last, end)-start+1;
      sublen = subend-substart;
      //cerr << "replacing bases " << substart << " through " << sublen << " with N" << endl;
      out.replace(substart,sublen,sublen,'N');
      //cerr << out << endl;
    }
  }
  
  //cerr << out << endl;
  return(out);
}

  
void TwoBit::printHeader(ostream & os)
{
  os << "Signature: " << hex << header.signature << endl;
  os << "Version:   " << dec << header.version << endl;
  os << "Seq Count: " << dec << header.sequenceCount << endl;
  os << "Reserved:  " << dec << header.reserved << endl << endl;
}

void TwoBit::printIndex(ostream & os)
{
  unsigned int size = index.size();
  for (unsigned int i=0; i<size; i++)
  {
    os << "Name:   " << index[i].name << endl
       << "Offset: " << index[i].offset << endl << endl;
  }
}




