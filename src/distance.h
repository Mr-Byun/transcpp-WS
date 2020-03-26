/*********************************************************************************
*                                                                                *
*     distance.h                                                                 *
*                                                                                *
*     Contains functions and interfaces for distance functions                   *
*                                                                                *
*********************************************************************************/

# ifndef DISTANCE_H
# define DISTANCE_H

#include "parameter.h"
#include "mode.h"

#include <map>
#include <boost/unordered_map.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <iostream>



/************************    Distance Functions   *******************************/

/*  All distance functions take a distance and return a value between 0 and 1 as
a function of that distance. They accept a vector of doubles as input, that way 
they can be pointed to by any of the interactions */

// returns 1 if less than or equal to max, else return 0
//double Uniform(double distance, double max);

// returns 1 up to a, then linear decay to b
//double Trapezoid(double distance, double a, double b);

// linear decay
//double Linear(double distance, double max);

// A linearly decaying sine wave
//double Sin(double distance, double period, double offset, double max);
  

/************************    Distance Class   ***********************************/

// interface and contain for distance functions and their parameters
class Distance
{
private:
  string name;
  string func_name;
  
  map<string, double_param_ptr> params;
  
  double max_distance;
  
  boost::function<double (double)> distFunc;
  
  mode_ptr mode;
  
public:
  Distance();
  Distance(ptree& pt);
  
  double getMaxDistance();
  double getDistFunc(double);
  string getName() { return(name); }
  void   getParameters(param_ptr_vector& p);
  void   getAllParameters(param_ptr_vector& p);
  map<string, double_param_ptr>& getParams() { return params; }
  
  void setParam(string, double);
  void setDistFunc(string funcname);
  void setMode(mode_ptr mode) { this->mode = mode; }
  
  void print(ostream& os);
  void read(ptree& pt);
  void write(ptree&);
};

typedef boost::shared_ptr<Distance> distance_ptr;

class DistanceContainer
{
private:
  map<string, distance_ptr> distances;
  mode_ptr mode;
  
  
public:
  DistanceContainer();
  DistanceContainer(ptree& pt);
  
  void getParameters(param_ptr_vector& p);
  void getAllParameters(param_ptr_vector& p);
  
  void add(string n, distance_ptr d) {distances[n] = d;} 
  void add(ptree& pt);
  void setMode(mode_ptr mode) { this->mode = mode; }
  
  void print(ostream& os);
  void write(ptree&);
  
  distance_ptr getDistance(string name);
  
  double distFunc(string name, double distance);
};
  
  


typedef boost::shared_ptr<DistanceContainer> distances_ptr;





#endif
  

