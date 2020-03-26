#include <Rcpp.h>
#include "organism.h"

#ifndef R_PWM_H
#define R_PWM_H

using namespace Rcpp;

class PWMPtr
{
private:
  PWM* pwm;
  
public:
  // constructors
  PWMPtr() : pwm(new PWM) {}
  PWMPtr(NumericMatrix mat, string type);
  PWMPtr(NumericMatrix mat, string type, double gc, double pseudo);
  PWMPtr(PWM* pwm) { this->pwm = pwm; } // cannot be exposed!
  
  // getters
  string get_source()  { return pwm->getSource();   }
  double get_gc()      { return pwm->getGC();       }
  double get_pseudo()  { return pwm->getPseudo();   }
  double get_max()     { return pwm->getMaxScore(); }
  int    get_length()  { return pwm->length();      }
  
  NumericMatrix get_pwm(string type);
  NumericMatrix get_pcm()  { return get_pwm("PCM");  }
  NumericMatrix get_pfm()  { return get_pwm("PFM");  }
  NumericMatrix get_pssm() { return get_pwm("PSSM"); }
  
  // setters
  void set_source(string source)  { pwm->setSource(source);   }
  void set_gc(double gc)          { pwm->setGC(gc);           }
  void set_pseudo(double pseudo)  { pwm->setPseudo(pseudo);   }
  
  void set_pwm(NumericMatrix mat, string type);

  // methods
  double pval2score(double pval)  { pwm->pval2score(pval);  }  // returns the threshold that would yeild a given p-value
  double score2pval(double score) { pwm->score2pval(score); } // returns the pvalue of a given score
  List   score(CharacterVector seq);
  
};






















#endif
