/*********************************************************************************
*                                                                                *
*     competition.h                                                              *
*                                                                                *
*     A very simple class holding the parameters for promoter competition        *
*                                                                                *
*********************************************************************************/


#ifndef COMP_H
#define COMP_H

#include "parameter.h"
#include "mode.h"
#include "utils.h"

class Competition
{
private:
  mode_ptr  mode;
  
  // basic sliding window parameters
  double_param_ptr window;      // the size of dna that can interact with a promoter
  double_param_ptr shift;       // the number of bases to shift to call a new window
  
  // parameters when using sum method
  double_param_ptr specificity; // and exponent that scales how much N affects promoter contact
  double_param_ptr threshold;   // the minimum N to have any probability of promoter contact
  double_param_ptr background;  // the relative probability that no peice of DNA interacts
  double_param_ptr S; // the free energy contribution of a single activator on DNA
  double_param_ptr interaction_strength;
  
  bool product; // whether T is proportional to the sum or product of Ns
  
  
public:
  // Constructors
  Competition();
  Competition(ptree& pt, mode_ptr mode);
  
  // Getters
  void getParameters(param_ptr_vector& p);
  void getAllParameters(param_ptr_vector& p);
  double getWindow()      { return window->getValue();      }
  double getShift()       { return shift->getValue();       }
  double getSpecificity() { return specificity->getValue(); }
  double getThreshold()   { return threshold->getValue();   }
  double getBackground()  { return background->getValue();  }
  double getS()           { return S->getValue();           }
  bool   getProduct()     { return product; }
  double getInteractionStrength() { return interaction_strength->getValue(); }
  
  // Setters
  void setWindow(double x)      { window->set(x);      }
  void setShift(double x)       { shift->set(x);       }
  void setSpecificity(double x) { specificity->set(x); }
  void setThreshold(double x)   { threshold->set(x);   }
  void setBackground(double x)  { background->set(x);  }
  void setS(double x)           { S->set(x);           }
  void setProduct(bool product) { this->product = product; }
  void setInteractionStrength(double x) { interaction_strength->set(x); }
  
  // I/O
  void set(ptree& pt, mode_ptr mode);
  void write(ptree& pt);
};

typedef boost::shared_ptr<Competition> competition_ptr;

#endif
