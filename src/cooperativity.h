/*********************************************************************************
*                                                                                *
*     cooperativity.h                                                            *
*                                                                                *
*     Holds information about cooporativity between tfs                          *
*                                                                                *
*********************************************************************************/

#ifndef COOP_H
#define COOP_H

#include "parameter.h"
#include "distance.h"
#include "mode.h"

#include <map>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>


/***********************    Cooporativity Class   *******************************/

class Cooperativity
{
private:
  // sources of information
  distances_ptr distances;
  mode_ptr      mode;
  
  // tfs involved
  string factor1;
  string factor2;
  
  // parameters
  double_param_ptr    Kcoop;
  distance_ptr dist;
  
  // Orientations
  bool HH;
  bool HT;
  bool TT;
  
public:
  // Constructors
  Cooperativity();
  Cooperativity(distances_ptr, ptree&);
  
  // Getters
  void   getParameters(param_ptr_vector& p);
  void   getAllParameters(param_ptr_vector& p);
  pair<string,string> getTFs();
  distance_ptr getDist() { return dist; }
  double distFunc(double d) { return dist->getDistFunc(d); }
  double getK() { return Kcoop->getValue(); }
  bool getHH() { return HH; }
  bool getHT() { return HT; }
  bool getTT() { return TT; }
  double_param_ptr getKcoopParam() { return Kcoop; }
  
  // Setters
  void setTFs(string a, string b) { factor1=a; factor2=b;}
  void setDist(distance_ptr d)    { dist = d;}
  void setK(double k)             { Kcoop->set(k);}
  void setMode(mode_ptr     mode) { this->mode = mode;}
    
  // I/O
  void read(ptree& pt);
  void write(ptree& pt);
  void print(ostream& os);
};
  
typedef boost::shared_ptr<Cooperativity> coop_ptr;
typedef vector<pair<string,coop_ptr> > coop_pairs;

/* A container which stores multiple Promoter objects and interfaces with
genes */
class CooperativityContainer
{
private:
  // sources of information
  distances_ptr distances;
  mode_ptr      mode;
  
  vector<coop_ptr> coops;
  
public:
  //Constructors
  CooperativityContainer();
  CooperativityContainer(distances_ptr, ptree&);
  
  // Getters
  coop_pairs   getCoops(string); // returns a pair of the tf this coops with and a pointer to the coop object
  void         getParameters(param_ptr_vector& p);
  void         getAllParameters(param_ptr_vector& p);
  vector<coop_ptr>& getAllCoops() { return coops; }
  
  // Setters
  void add(coop_ptr c) {coops.push_back(c);}
  void setDistances(distances_ptr d) {distances = d;}
  void setMode(mode_ptr mode) { this->mode = mode;}
  
  // I/O
  void read(ptree& pt, distances_ptr);
  void write(ptree& pt);
  void print(ostream& os);
};
  
typedef boost::shared_ptr<CooperativityContainer> coops_ptr;
  
  
  
#endif
