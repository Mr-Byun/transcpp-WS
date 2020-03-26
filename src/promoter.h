/*********************************************************************************
*                                                                                *
*     promoter.h                                                                 *
*                                                                                *
*     Contains methods for how N is interpreted to yield a rate                  *
*                                                                                *
*********************************************************************************/

#ifndef PROMOTER_H
#define PROMOTER_H

#include "parameter.h"
#include "mode.h"

#include <map>
#include <boost/unordered_map.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <iostream>


/************************   Functions   ****************************************/

// Sigmoid shaped diffusion limited Arrhenius rate law
//double Arrhenius(double M, double max, double theta, double Q);


/************************    Promoter Class   ***********************************/


class Promoter
{
private:
  string name;
  string func_name;
  
  map<string, double_param_ptr> params;
  mode_ptr mode;
  
  boost::function<double (double)> rateFunc;
  
public:
  // Constructors
  Promoter();
  Promoter(ptree& pt);
  
  // Getters
  string& getName() {return name;}
  void    getParameters(param_ptr_vector& p);
  void    getAllParameters(param_ptr_vector& p);
  double  getRate(double M) { return rateFunc(M); }
  map<string, double_param_ptr>& getParamMap() { return params; }
  
  // Setters
  void setMode(mode_ptr mode) { this->mode = mode; }
  
  // I/O
  void read(ptree& pt);
  void write(ptree& pt);
  void print(ostream& os);
};
  
typedef boost::shared_ptr<Promoter> promoter_ptr;
  
/* A container which stores multiple Promoter objects and interfaces with
genes */
class PromoterContainer
{
private:
  map<string, promoter_ptr> promoters;
  mode_ptr mode;
  
public:
  //Constructors
  PromoterContainer();
  PromoterContainer(ptree& pt);
  
  // Getters
  promoter_ptr getPromoter(string& name) {return promoters[name];}
  void         getParameters(param_ptr_vector& p);
  void         getAllParameters(param_ptr_vector& p);
  
  // Setters
  void setMode(mode_ptr mode) { this->mode = mode; }
  
  // I/O
  void read(ptree& pt);
  void write(ptree& pt);
  void print(ostream& os);
};
  
typedef boost::shared_ptr<PromoterContainer> promoters_ptr;
  
  
  
#endif
