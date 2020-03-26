/*********************************************************************************
*                                                                                *
*     chromatin.h                                                                *
*                                                                                *
*     Contains the accessibility state for each gene                             *
*                                                                                *
*********************************************************************************/

#ifndef CHROMATIN_H
#define CHROMATIN_H

#include <vector>
#include <map>
#include <boost/shared_ptr.hpp>

#include "gene.h"
#include "mode.h"

class Chromatin
{
private:
  genes_ptr genes;
  mode_ptr  mode;
  
  double_param_ptr kacc;
  
  map<Gene*, vector<double> > accessibility;
  
public:
  Chromatin();
  
  void read(ptree&,genes_ptr,mode_ptr);
  void write(ptree&);
  
  void setGenes(genes_ptr genes) { this->genes = genes; }
  void setGenes(mode_ptr  mode)  { this->mode  = mode;  }
  
  vector<double>& getAcc(Gene& gene) { return accessibility[&gene]; }
  double getKacc() {return kacc->getValue();}
  
  void getParameters(param_ptr_vector& p);
  void getAllParameters(param_ptr_vector& p);
  
};

typedef boost::shared_ptr<Chromatin> chromatin_ptr;

#endif