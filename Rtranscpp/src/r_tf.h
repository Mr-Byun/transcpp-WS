#include <Rcpp.h>
#include "r_pwm.h"
#include "organism.h"

#ifndef R_TF_H
#define R_TF_H

using namespace Rcpp;

class TFPtr
{
private:
  tf_ptr tf;
  
public:
  // constructors
  TFPtr() : tf(tf_ptr(new TF)) {}
  TFPtr(tf_ptr tf) { this->tf = tf; } // cannot be exposed!
  
  // getters
  vector<double> get_coefs()   { return tf->getCoefs();        }
  
  string get_name()            { return tf->getName();         }
  double get_lambda()          { return tf->getLambda();       }
  double get_threshold()       { return tf->getThreshold();    }
  int    get_binding_size()    { return tf->getBindingSize();  }

  PWMPtr get_pwm() { return PWMPtr(tf->getPWMPtr()); }
  
  // setters
  void set_name(string name)       { tf->setName(name);        }
  void set_source(string source)   { tf->setPWMSource(source); }
  void set_binding_size(int n)     { tf->setBindingSize(n);    }
  void set_coefs(vector<double> c) { tf->setCoefs(c);          }    
  void set_lambda(double l)        { tf->setLambda(l);         }
  void set_threshold(double t)     { tf->setThreshold(t);      }
  
};



























#endif
