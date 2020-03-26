#include "r_parameter.h"
#include "r_pwm.h"

using namespace Rcpp;

DoubleParameterPtr::DoubleParameterPtr(organism_ptr parent, int index)
{
  this->parent = parent;
  this->index  = index;
  param_ptr_vector& all_params = parent->getAllParameters();
  iparam_ptr iparam = all_params[index];
  param = boost::dynamic_pointer_cast<Parameter<double> >(iparam);
}

PWMParameterPtr::PWMParameterPtr(organism_ptr parent, int index)
{
  this->parent = parent;
  this->index  = index;
  param_ptr_vector& all_params = parent->getAllParameters();
  iparam_ptr iparam = all_params[index];
  param = boost::dynamic_pointer_cast<Parameter<vector<vector<double> > > >(iparam);
}

SeqParameterPtr::SeqParameterPtr(organism_ptr parent, int index)
{
  this->parent = parent;
  this->index  = index;
  param_ptr_vector& all_params = parent->getAllParameters();
  iparam_ptr iparam = all_params[index];
  param = boost::dynamic_pointer_cast<Parameter<Sequence> >(iparam);
}
  
NumericMatrix PWMParameterPtr::getTypeValue(string stype)
{
  vector<vector<double> >& pwm = param->getValue();
  
  int type;
  if (stype == string("PCM"))
    type = PCM;
  else if (stype == string("PFM"))
    type = PFM;
  else if (stype == string("PSSM"))
    type = PSSM;
  else 
    error("type must be one of PCM, PFM, PSSM"); 
  
  //vector<vector<double> > matrix = pwm.getPWM(type);
  int length = pwm.size();
  NumericMatrix out(length, 4);
  for (int i=0; i<length; i++)
  {
    for (int j=0; j<4; j++)
      out(i,j) = pwm[i][j];
  }
  return out;
}

void PWMParameterPtr::setTypeValue(NumericMatrix x, string stype)
{
  vector<vector<double> >& pwm = param->getValue();
  
  int type;
  if (stype == string("PCM"))
    type = PCM;
  else if (stype == string("PFM"))
    type = PFM;
  else if (stype == string("PSSM"))
    type = PSSM;
  else 
    error("type must be one of PCM, PFM, PSSM"); 
  
  int length = x.ncol();
  vector<vector<double> > matrix(length, vector<double>(4));

  for (int i=0; i<length; i++)
  {
    for (int j=0; j<4; j++)
      pwm[i][j] = x(i,j);
  }
  
  //pwm.setPWM(matrix, type);
  parent->move_all(index);
}

CharacterVector SeqParameterPtr::getValue()
{
  vector<char> out;
  vector<int>& seq = param->getValue().getSequence();
  int length = seq.size();
  out.resize(length);
  for (int i=0; i<length; i++)
  {
    switch (seq[i])
    {
    case 0:
      out[i] = 'A';
      break;
    case 1:
      out[i] = 'C';
      break;
    case 2:
      out[i] = 'G';
      break;
    case 3:
      out[i] = 'T';
      break;
    case 4:
      out[i] = 'N';
      break;
    }
  }
  CharacterVector output = wrap(out);
  return output;
}
  

void SeqParameterPtr::setValue(CharacterVector input)
{
  string out;
  CharacterVector::iterator i;
  for (i = input.begin(); i != input.end(); i++)
    out += *i;
  param->set(out);
  parent->move_all(index);
}

void DoubleParameterPtr::setValue(double input)
{
  param->set(input);
  parent->move_all(index);
}


CharacterVector DoubleParameterPtr::TFName()
{
  if (param->is_tf_param())
    return param->getTFName();
  else
    return NA_STRING;
}

CharacterVector PWMParameterPtr::TFName()
{
  if (param->is_tf_param())
    return param->getTFName();
  else
    return NA_STRING;
}

CharacterVector SeqParameterPtr::TFName()
{
  if (param->is_tf_param())
    return param->getTFName();
  else
    return NA_STRING;
}

RCPP_MODULE(mod_param)
{
  class_<DoubleParameterPtr>("DoubleParameter")
  .constructor()
  .property("value",  &DoubleParameterPtr::getValue, &DoubleParameterPtr::setValue)
  .property("tf",     &DoubleParameterPtr::TFName   )
  .property("name",   &DoubleParameterPtr::name     )
  .property("type",   &DoubleParameterPtr::type     )
  .property("anneal", &DoubleParameterPtr::getAnneal)
  
  .method("restore", &DoubleParameterPtr::restore)
  ;
  
  class_<PWMParameterPtr>("PWMParameter")
  .constructor()
  .property("value",  &PWMParameterPtr::getValue, &PWMParameterPtr::setValue)
  .property("tf",     &PWMParameterPtr::TFName   )
  .property("name",   &PWMParameterPtr::name     )
  .property("type",   &PWMParameterPtr::type     )
  .property("anneal", &PWMParameterPtr::getAnneal)
  
  .method("restore", &PWMParameterPtr::restore)
  ;
  
  class_<SeqParameterPtr>("SeqParameter")
  .constructor()
  .property("value",  &SeqParameterPtr::getValue, &SeqParameterPtr::setValue)
  .property("tf",     &SeqParameterPtr::TFName   )
  .property("name",   &SeqParameterPtr::name     )
  .property("type",   &SeqParameterPtr::type     )
  .property("anneal", &SeqParameterPtr::getAnneal)
  
  .method("restore", &SeqParameterPtr::restore)
  ;
}



