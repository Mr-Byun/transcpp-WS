#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "utils.h"

class Sequence 
{
private:
  vector<int> sequence;
public:
  Sequence();
  Sequence(const string& s);
  Sequence(const vector<int>& s);
  Sequence(const vector<char>& s);

  void setSequence(const string& s);
  void setSequence(const vector<int>& s);
  void setSequence(const vector<char>& s);

  vector<int>& getSequence();
  vector<int>& getSequence() const;
  string getSequenceString();
  int getLength() { return sequence.size(); }
  
  void   serialize(void *buf) const;
  void   deserialize(void const *buf);
  size_t getSize();
  
  void print(ostream& os);
  
};

#endif
