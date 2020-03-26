/*********************************************************************************
*                                                                                *
*     coeffects.h                                                                *
*                                                                                *
*     The structures for coactivation and corepression.                          *
*                                                                                *
*********************************************************************************/

#ifndef COEFFECTS_H
#define COEFFECTS_H

#include "parameter.h"
#include "distance.h"
#include "mode.h"

#include <map>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>


/***********************    Coeffect Class   *******************************/

class Coeffect
{
private:
  // sources of information
  distances_ptr distances;
  mode_ptr      mode;
  
  // tfs involved
  string actor;
  string target;
  int    coef_idx; // what coef of the target this causes
  
  // parameters
  double_param_ptr    efficiency;
  distance_ptr dist;
  
  // Orientations
  bool HH;
  bool HT;
  bool TT;
  
public:
  // Constructors
  Coeffect();
  Coeffect(distances_ptr, ptree&);
  
  // Getters
  void   getParameters(param_ptr_vector& p);
  void   getAllParameters(param_ptr_vector& p);
  distance_ptr getDist() {return dist; }
  
  string& getActor()      { return actor; }
  string& getTarget()     { return target; }
  double  getEfficiency() { return efficiency->getValue(); }
  int     getIdx()        { return coef_idx; }
  double_param_ptr getEfficiencyParam() { return efficiency; }
  
  double distFunc(double d) { return dist->getDistFunc(d); }
  double getMaxDistance()   { return dist->getMaxDistance(); }
  
  // Setters
  void setDist(distance_ptr dist) { this->dist = dist;}
  void setMode(mode_ptr     mode) { this->mode = mode;}
    
  // I/O
  void read(ptree& pt);
  void write(ptree& pt);
  void print(ostream& os);
};
  
typedef boost::shared_ptr<Coeffect> coeffect_ptr;

typedef vector< pair<string,coeffect_ptr> > coeffect_pairs;

/* A container which stores multiple Promoter objects and interfaces with
genes */
class CoeffectContainer
{
private:
  // sources of information
  distances_ptr distances;
  mode_ptr      mode;
  
  vector<coeffect_ptr> coeffects;
  
public:
  //Constructors
  CoeffectContainer();
  CoeffectContainer(distances_ptr, ptree&);
  
  // Getters
  coeffect_pairs getTargets(string tfname);
  void           getParameters(param_ptr_vector& p);
  void           getAllParameters(param_ptr_vector& p);
  vector<coeffect_ptr>& getAllCoeffects() { return coeffects; }
  
  // Setters
  void add(coeffect_ptr c) {coeffects.push_back(c);}
  void setDistances(distances_ptr d) {distances = d;}
  void setMode(mode_ptr mode) { this->mode = mode;}
  
  // I/O
  void read(ptree& pt, distances_ptr);
  void write(ptree& pt);
  void print(ostream& os);
};
  
typedef boost::shared_ptr<CoeffectContainer> coeffects_ptr;
  
  
  
#endif
