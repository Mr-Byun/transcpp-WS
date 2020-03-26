/*********************************************************************************
*                                                                                *
*     scalefactor.h                                                             *
*                                                                                *
*     Contains a simple structure for scale factors                              *
*                                                                                *
*********************************************************************************/

#ifndef SCALEFACTOR_H
#define SCALEFACTOR_H

#include "parameter.h"
#include "mode.h"

#include <map>
#include <boost/unordered_map.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>

using namespace std;

class ScaleFactor
{
private:
  /* I use an affine ( ax + b ) transformation for scaling data. I generally
  assume that b=0, but some users set offset lower or are more or less 
  conservative in background removal */
  double_param_ptr A;
  double_param_ptr B;
  
  string name;
  
  mode_ptr mode;
  
public:
  // Constructors
  ScaleFactor();
  ScaleFactor(ptree& pt);
  
  // Getters
  void getParameters(param_ptr_vector& p);
  void getAllParameters(param_ptr_vector& p);
  double_param_ptr getA()  { return A;    }
  double_param_ptr getB()  { return B;    }
  string& getName() { return name; }
  
  // Setters
  void setMode(mode_ptr     mode) { this->mode = mode;}
  
  // Method
  double scale(double x);
  double unscale(double x);
  
  // I/O
  void read(ptree& pt);
  void write(ptree& pt);
};

typedef boost::shared_ptr<ScaleFactor> scale_factor_ptr;

class ScaleFactorContainer
{
private:
  map<string, scale_factor_ptr> scales;
  vector<scale_factor_ptr> scales_vector; // simply easier to loop over
  mode_ptr mode;
  
public:
  // Constructors
  ScaleFactorContainer();
  ScaleFactorContainer(ptree& pt);
  
  // Getters
  void             getParameters(param_ptr_vector& p);
  void             getAllParameters(param_ptr_vector& p);
  scale_factor_ptr getScaleFactor(string name);
  
  // Setters
  void setMode(mode_ptr     mode) { this->mode = mode;}
  
  // I/O
  void read(ptree& pt);
  void write(ptree& pt);
};

typedef boost::shared_ptr<ScaleFactorContainer> scale_factors_ptr;

#endif
