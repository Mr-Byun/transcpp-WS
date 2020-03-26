/*********************************************************************************
*                                                                                *
*     parameter.h                                                                *
*                                                                                *
*     Contains class description for parameters                                  *
*                                                                                *
*********************************************************************************/


#ifndef PARAMETER_H
#define PARAMETER_H

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>

#include "sequence.h"

using namespace std;
using boost::property_tree::ptree;


/* 
Ok, this class got a lot more sophisticated! Parameters are now templated
and abstracted. The basic functions are handled by an interface class that serves
as a base class to all the templated classes. This way, we can hold a vector of
pointers to the individual parameters using the interface class and tweak can
be called without knowing the actual type of the parameter.

To add new parameter types the user must override the functions that operate on
the value itself, such as the get and set value functions, checkLimits, and tweak
*/

// a boost uniform double distribution for generating variates. More robust than
// basic c random functions.
class UnitDblDist
{
private:
  boost::minstd_rand baseGen;
  boost::uniform_real<> uniDblUnit;
  boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uniDblGen;
  
public:
  UnitDblDist():
    baseGen(0),uniDblUnit(0,1),uniDblGen(baseGen, uniDblUnit) {}
  
  double draw() { return uniDblGen(); }
};
    

class ParameterInterface
{
protected:
  ptree* node;            // points to the node used to create this (useful in scrambling)

  string param_name;
  string type;
  string tf_name;
  string move_func;
  string restore_func;

  bool out_of_bounds;     // is it outside of the limits?
  bool anneal;            // are we ever tweaking it?
  
  bool tf_name_set;
  
  //unsigned int seed;      // seed used for random number generation
  
  virtual bool checkLimits() = 0;
  void setTypeName();
  
  // random number generation
  UnitDblDist dist;
  
public:
  
  virtual ~ParameterInterface() {}
  
  // Getters
  bool isOutOfBounds() { return out_of_bounds; }
  bool isAnnealed()    { return anneal;        }
  bool is_tf_param()   { return tf_name_set;   }

  const string& getParamName() { return param_name;   }
  const string& getTFName()    { return tf_name;      }
  const string& getMove()      { return move_func;    }
  const string& getRestore()   { return restore_func; }  
  const string& getType()      { return type;         }
  
  // Setters
  void setParamName(const string& s) { this->param_name   = s;   }
  void setTFName(const string& s)    { this->tf_name      = s; tf_name_set = true; }
  void setAnnealed(bool b)           { this->anneal       = b;   }
  void setMove(string s)             { this->move_func    = s;   }
  void setRestore(string s)          { this->restore_func = s;   }
  void setType(string type)          { this->type = type;        }
  void setNode(ptree* node)          { this->node = node;        }
  
  // Methods
  virtual void scramble(double) = 0;
  virtual void tweak(double)    = 0;
  virtual void restore()        = 0;
  
  // Serialize
  virtual void   serialize(void *buf) const   = 0;
  virtual void   deserialize(void const *buf) = 0;
  virtual size_t getSize()                    = 0;
  
  // I/O
  virtual void write(ptree& pt, int precision) const   = 0;
  virtual void print(ostream& os)                      = 0;
  virtual void printHeader(ostream& os)                = 0;
};
  
template< typename T> 
class Parameter : public ParameterInterface
{
private:
  T value;           // the value of the parameter
  T previous_value;  // the last value the parameter had
  
  double lim_low;
  double lim_high;
 
  bool checkLimits();
  void setTypeName();
  
public:
  // Constructor 
  Parameter();
  Parameter(string s) { param_name = s; }
  Parameter(ptree& pt);
  Parameter(string, ptree& pt);
  Parameter(ptree& pt, string name, string move);
  Parameter(string name, string move);
  
  void read(ptree& pt);
  
  // Getters
  T& getValue();
  T getPrevious() const;
  
  double getLimHigh() const;
  double getLimLow() const;
  
  // Setters
  void set(T);
  void setLimits(double low, double high);
  
  // Methods
  void scramble(double);
  void tweak(double delta);
  void restore();
  
  // Serialize
  void   serialize(void *buf) const;
  void   deserialize(void const *buf);
  size_t getSize();                    
  
  // Output
  void write(ptree& pt, int precision) const;
  void print(ostream& os);
  void printHeader(ostream& os);
};

typedef boost::shared_ptr<Parameter<double>   > double_param_ptr;
typedef boost::shared_ptr<Parameter<int>      > int_param_ptr;
typedef boost::shared_ptr<Parameter<Sequence> > seq_param_ptr;
typedef boost::shared_ptr<Parameter<vector<vector<double> > > > pwm_param_ptr;
typedef boost::shared_ptr<ParameterInterface>   iparam_ptr;

typedef vector<iparam_ptr> param_ptr_vector;


#endif
