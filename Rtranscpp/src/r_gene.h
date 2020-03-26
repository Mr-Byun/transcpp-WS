#include <Rcpp.h>
#include "organism.h"

#ifndef R_GENE_H
#define R_GENE_H

class GenePtr
{
private:
  gene_ptr gene;
  
public:
  // constructors
  GenePtr() : gene(gene_ptr(new Gene)) {}
  GenePtr(gene_ptr gene) { this->gene = gene; } // cannot be exposed!
  
  // getters
  Rcpp::CharacterVector get_sequence(); 
  int          get_right_bound() { return gene->getRightBound();     }
  int          get_left_bound()  { return gene->getLeftBound();      }
  int          get_length()      { return gene->length();            }
  string       get_name()        { return gene->getName();           }
  gene_ptr     get_ptr()         { return gene;                      }
  
  // setters
  void set_sequence(Rcpp::CharacterVector input);
  
  
};



























#endif
