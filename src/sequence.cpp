#include "sequence.h"


Sequence::Sequence() {}

Sequence::Sequence(const string& s)       { sequence = string2int(s); }

Sequence::Sequence(const vector<int>& s)  { sequence = s; }

Sequence::Sequence(const vector<char>& s) { sequence = char2int(s); }


void Sequence::setSequence(const string& s)      { sequence = string2int(s); }
void Sequence::setSequence(const vector<int>& s) { sequence = s; }

vector<int>&  Sequence::getSequence() {return sequence; }

string Sequence::getSequenceString() {return int2string(sequence); }


void Sequence::serialize(void *buf) const
{
  int n = sequence.size();
  
  int* dest = static_cast<int *>(buf);
  dest[0] = n;
  for (int i=0; i<n; i++)
    dest[i+1] = sequence[i];
}

void Sequence::deserialize(void const *buf)
{
  int const * from = static_cast<int const *>(buf);
  int n = from[0];
  sequence.resize(n);
  for (int i=0; i<n; i++)
    sequence[i] = from[i+1];
}

size_t Sequence::getSize()
{
  int n = sequence.size() + 1;
  return n*sizeof(int);
}

void Sequence::print(ostream& os)
{
  string seq = int2string(sequence);
  os << seq;
}
  
