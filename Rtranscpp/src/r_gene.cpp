#include "r_gene.h"

using namespace Rcpp;

CharacterVector GenePtr::get_sequence()
{
  CharacterVector output = wrap(gene->getSequenceChars());
  return output;
}
  

void GenePtr::set_sequence(CharacterVector input)
{
  seq_param_ptr sequence = gene->getSequenceParam();
  string out;
  CharacterVector::iterator i;
  for (i = input.begin(); i != input.end(); i++)
    out += *i;
  sequence->set(out);
}
  
  
  
  
  
RCPP_MODULE(mod_gene)
{
  class_<GenePtr>("Gene")
  
  .constructor()
  
  .property("sequence",    &GenePtr::get_sequence    , &GenePtr::set_sequence) 
  .property("right_bound", &GenePtr::get_right_bound  )
  .property("left_bound",  &GenePtr::get_left_bound   ) 
  .property("length",      &GenePtr::get_length       )     
  .property("name",        &GenePtr::get_name         )       
  
  ;
}
