/*********************************************************************************
*                                                                                *
*     TF.h                                                                       *
*                                                                                *
*     Contains class definition for Transcription Factors                        *
*                                                                                *
*********************************************************************************/

#ifndef TF_H
#define TF_H

#include "parameter.h"
#include "pwm.h"
#include "cooperativity.h"
#include "coeffects.h"
#include "mode.h"
#include <vector>
#include <string>
#include <cstdlib>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>



using namespace std;
using boost::property_tree::ptree;

class TF : boost::noncopyable
{
private: 
  mode_ptr  mode;
  
  int    index;     // which index is this TF within an organism?
                    
  PWM       energy; // binding energy model
  string    tfname; 
  int       bsize;  // the footprint size of the TF
  int       offset; // the offset of the motif within a footprint (NOT IMPLEMENTED YET)
                       
  double_param_ptr kmax;       // true binding affinity confounded with relative protein levels
  double_param_ptr threshold;  // a somewhat arbitrary cutoff to use for calling a binding site
  double_param_ptr lambda;     // how much a mutation away from consensus effects ddG
  
  double Kns; // the nonspecific binding energy
  
  vector<double_param_ptr> coefs;
  
  // index is -log10(pval)*precision
  vector<double> s2p;

  vector< pair<TF*, coop_ptr> >     coops;
  vector< pair<TF*, coeffect_ptr> > coeffects;
  
  //void subscore(const vector<int> & s, double* out); // for scoring sequence of pwmlen
  void readPWM(ptree& pt);
  
public:
  // constructors
  TF();           // empty constructor
  TF(ptree& pt, mode_ptr);  // create from property tree
  
  // getters
  const string&  getName() const;
  const string&  getPWMSource() {return energy.getSource();}
  int            getPWMLength() {return energy.length();   }
  double         getThreshold();
  double         getKns() {return Kns;}
  int            getBindingSize() const;
  double         getKmax();
  double&        getCoef();
  double         getModifiedCoef();
  double         getLambda();
  double         getMaxScore();
  void           getParameters(param_ptr_vector& p); // pushes parameters onto argument
  void           getAllParameters(param_ptr_vector& p); // pushes parameters onto argument
  bool           ncoops() { return coops.size(); } // used to check if this ever cooperates
  bool           checkCoops(TF*, char, char); // check if this cooperates with tf
  coop_ptr       getCoop(TF*);
  bool           neverActivates();
  bool           neverQuenches();
  vector<double> getCoefs();
  int            getNumModes() { return coefs.size(); }
  int            getIndex()    { return index; }
  
  vector<double_param_ptr>& getCoefPtrs() {return coefs; }
  double_param_ptr getThresholdParam() { return threshold; }
  double_param_ptr getLambdaParam()    { return lambda; }
  double_param_ptr getKmaxParam()      { return kmax; }
  pwm_param_ptr    getPWMParam()       { return energy.getPWMParam(); }
  PWM*             getPWMPtr()         { return &energy; }
  PWM&             getPWM()            { return energy; }
  
  vector< pair<TF*, coeffect_ptr> >& getTargets() { return coeffects; }   
  
  // setters
  void setPWMSource(string source) { energy.setSource(source); }
  void setNscore();
  void set(ptree& pt); // set entire TF
  void setName(string n);
  void setPWM(vector<vector<double> >& t, int type);
  void setPWM(vector<vector<double> >& t, int type, double gc, double pseudo);
  void setKmax(double k);
  void setLambda(double l);
  void setThreshold(double t);
  void setBindingSize(int b);
  void setCoops(vector< pair<TF*,coop_ptr> > p) { coops = p; }
  void setCoeffects(vector< pair<TF*, coeffect_ptr> > p) { coeffects = p; }
  void setIndex(int index) { this->index = index; }
  void setCoefs(vector<double>);
  void setMaxScore(double mscore) { energy.setMaxScore(mscore); } // WSB
  
  // methods
  TFscore score(const string & s);     // score a string with tf
  TFscore score(const vector<int>& s); // score an int vector (faster) 
  void    score(const vector<int>& s, TFscore &t); // pass by reference (fastest) , 
  
  // I/O
  void print(ostream& os);
  void write(ostream& os); // write XML to output
  void write(ptree& tfsnode); // write XML to property tree

};

typedef boost::shared_ptr<TF> tf_ptr;
typedef vector<tf_ptr> tf_ptr_vector;

class TFContainer
{
private:
  tf_ptr_vector tfs;
  mode_ptr      mode;

public:
  // Constructor
  TFContainer();
  TFContainer(ptree& pt, mode_ptr);
  
  // Getters
  TF& getTF(const string& n);
  TF& getTF(int index);
  tf_ptr getTFptr(int index);
  tf_ptr getTFptr(const string& n);
  int size() const;
  void getParameters(param_ptr_vector& p);
  void getAllParameters(param_ptr_vector& p);
  
  // Setters
  void add(ptree& pt, mode_ptr);
  void add(tf_ptr t);
  void setCoops(coops_ptr);
  void setCoeffects(coeffects_ptr);
  
  // Output
  void write(ostream& os) const;
  void write(ptree& pt) const;
};

typedef boost::shared_ptr<TFContainer> tfs_ptr;

#endif
