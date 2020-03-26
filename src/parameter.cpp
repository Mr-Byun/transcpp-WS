/*********************************************************************************
*                                                                                *
*     parameter.cpp                                                              *
*                                                                                *
*     Contains class description for parameters                                  *
*                                                                                *
*********************************************************************************/


#include "parameter.h"
#include "utils.h"

#include <cmath>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

using boost::property_tree::ptree;

# define foreach_ BOOST_FOREACH

/*    Constructors    */

template< typename T> 
Parameter<T>::Parameter() 
{
  anneal        = false;
  tf_name_set   = false;
  param_name    = "not set";
  move_func     = "ResetAll";
  out_of_bounds = false;
  setTypeName();
}

template< typename T> 
Parameter<T>::Parameter(string name, string move) 
{
  anneal        = false;
  tf_name_set   = false;
  param_name    = name;
  move_func     = move;
  out_of_bounds = false;
  setTypeName();
}

template< typename T> 
Parameter<T>::Parameter(ptree& pt)
{
  node = &pt;
  setTypeName();
  read(pt);
  out_of_bounds = false;
}

template< typename T> 
Parameter<T>::Parameter(string name, ptree& pt)
{
  node       = &pt;
  param_name = name;
  setTypeName();
  read(pt);
}

// this is the prefered constructor!!!!
template< typename T> 
Parameter<T>::Parameter(ptree& pt, string name, string move)
{
  node = &pt;
  setTypeName();
  read(pt);
  out_of_bounds = false;
  param_name = name;
  move_func  = move;
}

template< typename T >
void Parameter<T>::setTypeName() { warning("typenames must be set to export parameters to R or Matlab"); }

template< >
void Parameter<int>::setTypeName()      { type = "int";      }

template< >
void Parameter<double>::setTypeName()   { type = "double";   }

template< >
void Parameter<Sequence>::setTypeName() { type = "Sequence"; }

template< >
void Parameter<vector<vector<double> > >::setTypeName()      { type = "PWM";      }


/*    Getters   */

template< typename T> 
T& Parameter<T>::getValue() { return value; }

template< typename T> 
T Parameter<T>::getPrevious() const  { return previous_value; }

template< typename T> 
double Parameter<T>::getLimHigh() const { return lim_high; }

template< typename T> 
double Parameter<T>::getLimLow() const { return lim_low; }


/*    Setters   */

template< typename T> 
void Parameter<T>::set(T v)
{
  previous_value = value;
  value = v;
}

template< typename T> 
void Parameter<T>::setLimits(double low, double high)
{
  lim_high = high;
  lim_low  = low;
}

/*    Methods   */

template< typename T> 
bool Parameter<T>::checkLimits()
{
  if (value > lim_high || value < lim_low)
  {
    out_of_bounds = true;
    return true;
  }
  else
  {
    out_of_bounds = false;
    return false;
  }
}

template<> 
bool Parameter<Sequence>::checkLimits()
{
  out_of_bounds = false;
  return false;
}

template<> 
bool Parameter<vector<vector<double> > >::checkLimits()
{
  return false;
}



template< typename T> 
void Parameter<T>::tweak(double delta)
{
  previous_value = value;
  value += delta;
  //cerr << "value moved from " << previous_value << " to " << value << endl;
  checkLimits();
}

/*  Naively, to tweak sequence I will simply pick some number of bases, based
on an integer rounding of delta, then mutate them to a random base. I will
use delta as a seed for random  */

template<> 
void Parameter<Sequence>::tweak(double delta)
{
  previous_value = value;
  
  if (delta < 0) delta = -delta;
  
  int nedits = ceil(delta); // the number of bases to mutate
  vector<int>& seq = value.getSequence(); // the seq to tweak
  int length = seq.size(); // the last index in the sequence
  
  if (nedits > length)
    nedits = length;
  if (nedits < 1)
    nedits = 1;
  
  for (int i=0; i<nedits; i++)
  {
    // output = min + (rand_r(seed) % (int)(max - min + 1))
    int pos  = dist.draw() * length;
    int base = dist.draw() * 4;
    seq[pos] = base;
  }
}

template<> 
void Parameter<vector<vector<double> > >::tweak(double delta)
{
  // we should only ever tweak from a move that is in bounds
  out_of_bounds = false;
  previous_value = value;
  
  int length = value.size();

  //for (int i=0; i<length; i++)
  //  cerr << value.getConsensus()[i];
  //cerr << endl;
  
  int pos1 = dist.draw() * length;
  int pos2 = dist.draw() * 4;
  
  value[pos1][pos2] += delta;
}

template< typename T> 
void Parameter<T>::scramble(double rand_uniform)
{
  value = (lim_high - lim_low)*rand_uniform + lim_low;
  stringstream tmp;
  tmp << setprecision(5);
  tmp.str("");
  tmp << value;
  node->put("<xmlattr>.value", tmp.str());
  if (checkLimits())
  {
    stringstream err;
    err << "ERROR: scrambled variable out of bounds" << endl;
    error(err.str());
  }
}

template<> 
void Parameter<Sequence>::scramble(double rand_uniform)
{
  vector<int>& seq = value.getSequence();
  int length = seq.size();
  
  //seed = (unsigned int) rand_uniform * 100000;
  
  for (int i=0; i<length; i++)
    seq[i] = dist.draw() * 4;
  
  string char_seq = int2string(seq);
  node->put("<xmlattr>.sequence", char_seq);
}

// I am making scrambled pwms according to 
template<> 
void Parameter<vector<vector<double> > >::scramble(double rand_uniform)
{
  stringstream tmp;
  //seed = (unsigned int) (rand_uniform * 100000);
  node->put("<xmlattr>.type","PSSM");
  
  int pwmpos = 0;
  foreach_(ptree::value_type& v, *node)
  {
    if (v.first == "position")
    {
      tmp.str("");
  
      for (int i=0; i<4; i++)
        value[pwmpos][i] = -10 + dist.draw() * 20;
      
      tmp << setw(10) << value[pwmpos][0] << ";"
          << setw(10) << value[pwmpos][1] << ";"
          << setw(10) << value[pwmpos][2] << ";"
          << setw(10) << value[pwmpos][3];
      v.second.put("",tmp.str());
      
      pwmpos++;
    }
  }
}
  


/*  Serialize and Deserialize functions  */

template< typename T> 
void Parameter<T>::serialize(void *buf) const
{
  T * dest = static_cast<T *>(buf);
  *dest = value;
}

template< typename T> 
void Parameter<T>::deserialize(void const *buf) 
{
  T const * from = static_cast<T const *>(buf);
  value = *from;
}

// we need to store sequence length as well, as this may change!
template<> 
void Parameter<Sequence>::serialize(void *buf) const
{
  value.serialize(buf); 
  //error("serialize not implemented for parameter of type Sequence");
}

template<> 
void Parameter<Sequence>::deserialize(void const *buf)
{

  value.deserialize(buf);
  //error("deserialize not implemented for parameter of type Sequence");
}


template<> 
void Parameter<vector<vector<double> > >::serialize(void *buf) const
{
  //value.serialize(buf); 
  error("serialize not implemented for parameter of type PWM");
}

template<> 
void Parameter<vector<vector<double> > >::deserialize(void const *buf)
{
  //value.deserialize(buf);
  error("deserialize not implemented for parameter of type PWM");
}


template< typename T> 
size_t Parameter<T>::getSize()
{
  return sizeof(T);
}

template<> 
size_t Parameter<Sequence>::getSize()
{
  return value.getSize();
}

template<> 
size_t Parameter<vector<vector<double> > >::getSize()
{
  //return value.getSize();
  error("getSize not implemented for parameter of type PWM");
}


template< typename T> 
void Parameter<T>::restore()
{
  value   = previous_value;
}


/*    I/O   */


template< typename T> 
void Parameter<T>::read(ptree& pt)
{
  node        = &pt;
  tf_name_set = false;
  
  // the move function is "ResetAll" my default. If it is defined here, override it
  move_func    = pt.get<string>("<xmlattr>.move", move_func);
  
  /* // need to let parameter get mode pointer 
  if (mode->getVerbose() >= 1 && (move_func == string("ResetAll") || restore_func == string("ResetAll")))
    cerr << "WARNING: A move or restore function was set to ResetAll. This can be very slow" << endl;
  */
   
  value        = pt.get<T>("<xmlattr>.value");
  anneal       = pt.get<bool>(  "<xmlattr>.anneal");
               
  lim_low      = pt.get<double>("<xmlattr>.lim_low");
  lim_high     = pt.get<double>("<xmlattr>.lim_high");
  
  previous_value = value;
  
  out_of_bounds = checkLimits();
}

template<> 
void Parameter<Sequence>::read(ptree& pt)
{
  stringstream err;
  err << "Reading of sequence is handled by gene.cpp!" << endl;
  error(err.str());
}

template<> 
void Parameter<vector<vector<double> > >::read(ptree& pt)
{
  stringstream err;
  err << "Reading of pwm is handled by TF.cpp!" << endl;
  error(err.str());
}


template< typename T> 
void Parameter<T>::write(ptree& pt, int precision) const
{
  stringstream tmp;
  tmp << setprecision(precision);
  tmp.str("");
  tmp << value;
  pt.put("<xmlattr>.value", tmp.str());
  tmp.str("");
  tmp << lim_low;
  pt.put("<xmlattr>.lim_low",  tmp.str());
  tmp.str("");
  tmp << lim_high;
  pt.put("<xmlattr>.lim_high", tmp.str());
  pt.put("<xmlattr>.anneal", anneal);
  pt.put("<xmlattr>.move", move_func);
}

template<> 
void Parameter<Sequence>::write(ptree& pt, int precision) const
{
  stringstream err;
  err << "ERROR: writing sequence not implemented yet" << endl;
  error(err.str());
}

template<> 
void Parameter<vector<vector<double> > >::write(ptree& pt, int precision) const
{
  stringstream err;
  err << "ERROR: writing PWM not implemented yet" << endl;
  error(err.str());
}

template< typename T> 
void Parameter<T>::printHeader(ostream& os)
{
  os << setprecision(4)
     << setw(10) << "TF"
     << setw(18) << "type"
     << setw(12) << "value"
     << setw(6)  << "tweak"
     << setw(10) << "lim_low"
     << setw(10) << "lim_high"
     << endl;
}

template< typename T> 
void Parameter<T>::print(ostream& os)
{
  os << setprecision(4)
     << setw(10) << tf_name
     << setw(18) << param_name
     << setw(12) << value
     << setw(6)  << anneal
     << setw(12) << lim_low
     << setw(12) << lim_high
     << endl;
}

template<> 
void Parameter<Sequence>::print(ostream& os)
{
  os << "    ";
  value.print(os);
}

template<> 
void Parameter<vector<vector<double> > >::print(ostream& os)
{
  //value.print(os,2);
  os << setprecision(4)
     << setw(10) << tf_name
     << setw(18) << param_name
     << setw(12) << "[4x" << value.size() << "]"
     << setw(6)  << anneal
     << setw(12) << "NA"
     << setw(12) << "NA"
     << endl;
}

template class Parameter<double>;
template class Parameter<int>;
template class Parameter<Sequence>;
template class Parameter<vector<vector<double> > >;


  
