#include <Rcpp.h>
#include "r_pwm.h"
#include "organism.h"

/* screw templates! They dont play well with Rcpp */

#ifndef R_PARAMETER_H
#define R_PARAMETER_H

using namespace Rcpp;

typedef boost::shared_ptr<Organism> organism_ptr;

class DoubleParameterPtr
{
private:
  int              index;
  organism_ptr     parent;
  double_param_ptr param;
  
public:
  DoubleParameterPtr() : param(double_param_ptr(new Parameter<double>)) {}
  DoubleParameterPtr(organism_ptr parent, int index);
  
  double getValue() { return param->getValue(); }
  void   setValue(double);
  
  void restore() { param->restore(); parent->restore_all(index); }
  
  bool getAnneal() { return param->isAnnealed(); }
  
  CharacterVector TFName();
  CharacterVector name() {return param->getParamName(); } 
  CharacterVector type() {return param->getType(); } 
};

class PWMParameterPtr
{
private:
  int              index;
  organism_ptr     parent;
  pwm_param_ptr    param;
  
public:
  PWMParameterPtr() : param(pwm_param_ptr(new Parameter<vector<vector<double> > >)) {}
  PWMParameterPtr(organism_ptr parent, int index);
  
  NumericMatrix getTypeValue(string type);
  void setTypeValue(NumericMatrix, string type);
  
  NumericMatrix getValue() { return getTypeValue("PSSM"); }
  void setValue(NumericMatrix x) { setTypeValue(x, "PSSM"); }
  
  void restore() { param->restore(); parent->restore_all(index); }
  
  bool getAnneal() { return param->isAnnealed(); }
  
  CharacterVector TFName();
  CharacterVector name() {return param->getParamName(); } 
  CharacterVector type() {return param->getType(); } 
};

class SeqParameterPtr
{
private:
  int              index;
  organism_ptr     parent;
  seq_param_ptr    param;
  
public:
  SeqParameterPtr() : param(seq_param_ptr(new Parameter<Sequence>)) {}
  SeqParameterPtr(organism_ptr parent, int index);
  
  CharacterVector getValue(); 
  void setValue(CharacterVector x);
  
  void restore() { param->restore(); parent->restore_all(index); }
  
  bool getAnneal() { return param->isAnnealed(); }
  
  CharacterVector TFName();
  CharacterVector name() {return param->getParamName(); } 
  CharacterVector type() {return param->getType(); } 
};
  



#endif


