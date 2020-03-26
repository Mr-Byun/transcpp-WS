/*********************************************************************************
*                                                                                *
*     distance.cpp                                                               *
*                                                                                *
*     Contains functions and interfaces for distance functions                   *
*                                                                                *
*********************************************************************************/

#include "distance.h"
#include "utils.h"

#include <math.h>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/ref.hpp>
#include <boost/unordered_map.hpp>

using namespace std;
using boost::property_tree::ptree;

#include <boost/foreach.hpp>
#include <limits>

# define foreach_ BOOST_FOREACH

/************************    Distance Functions   *******************************/

/*  All distance functions take a distance and return a value between 0 and 1 as
a function of that distance. They accept a vector of doubles as input, that way 
they can be pointed to by any of the interactions */

double Uniform(double distance, double max)
{
  if (distance < 0)
    distance = -distance;
  if (distance <= max)
    return 1;
  else
    return 0;
}

// a is the uniform distance, b is the linear distance after a
double Trapezoid(double distance, double a, double b)
{
  if (distance < 0)
    distance = -distance; 
  
  if (distance <= a)
    return 1;
  else if (distance < (a+b) )
  {
    double x   = distance      - a;
    return (1 - x/b);
  } else
    return 0;
}

double Linear(double distance, double max)
{
  if (distance < 0)
    distance = -distance; 
  
  if (distance < max)
    return (1 - distance/max);
  else
    return 0;
}

// a decaying sin wave
double Sine(double distance, double period, double offset, double max)
{
  const double pi = 3.1415926535897;
  double a      = 2*pi/period;
  
  if (distance < 0)
    distance = -distance; 
  
  if (distance < max)
    return (sin(a*distance+offset)+1)*(1-distance/max)/2;
  else
    return 0;
}


/************************    Distance Class   ***********************************/
    


/*    Constructors    */

Distance::Distance() {}

Distance::Distance(ptree& pt) { read(pt); }


/*    Getters   */

void Distance::getParameters(param_ptr_vector& p)
{
  typedef map<string, double_param_ptr>::iterator i_type;
  for(i_type i=params.begin(); i != params.end(); i++)
  {
    double_param_ptr& param = i->second;
    if (param->isAnnealed())
      p.push_back(param);
  }
}

void Distance::getAllParameters(param_ptr_vector& p)
{
  typedef map<string, double_param_ptr>::iterator i_type;
  for(i_type i=params.begin(); i != params.end(); i++)
  {
    double_param_ptr& param = i->second;
    p.push_back(param);
  }
}


/*    Setters   */

void Distance::setParam(string s, double v)
{
  if (params.find(s) == params.end())
  {
    double_param_ptr cur_param(new Parameter<double>());
    cur_param->set(v);
    params[s] = cur_param;
  } 
  else
    params[s]->set(v);
}

void Distance::setDistFunc(string funcname)
{
  if      (funcname == "Linear")
    distFunc = bind(&Linear, _1, boost::ref(params["Max"]->getValue()));
  
  else if (funcname == "Trapezoid")
    distFunc = bind(&Trapezoid, _1, 
      boost::ref(params["A"]->getValue()), 
      boost::ref(params["B"]->getValue()));
    
  else if (funcname == "Uniform")
    distFunc = bind(&Uniform, _1, boost::ref(params["Max"]->getValue()));
  
  else if (funcname == "Sine")
    distFunc = bind(&Sine, _1, 
      boost::ref(params["Period"]->getValue()), 
      boost::ref(params["Offset"]->getValue()), 
      boost::ref(params["Max"]->getValue()));
}
      
void Distance::read(ptree& pt)
{
  name      = pt.get<string>("<xmlattr>.name");
  func_name = pt.get<string>("<xmlattr>.distfunc");
  
  if      (func_name == "Linear")
  {
    double_param_ptr maxparam(new Parameter<double>(string(name+" Max"), pt.get_child("Max")));
    params["Max"] = maxparam;
    max_distance = params["Max"]->getLimHigh();
    distFunc = bind(&Linear, _1, boost::ref(params["Max"]->getValue()));

  } 
  else if (func_name == "Trapezoid")
  {
    double_param_ptr Aparam(new Parameter<double>(string(name+" A"), pt.get_child("A")));
    double_param_ptr Bparam(new Parameter<double>(string(name+" B"), pt.get_child("B")));
    params["A"] = Aparam;
    params["B"] = Bparam;
    max_distance = params["A"]->getLimHigh()+params["B"]->getLimHigh();
    distFunc = bind(&Trapezoid, _1, 
      boost::ref(params["A"]->getValue()), 
      boost::ref(params["B"]->getValue()));

  }
  else if (func_name == "Uniform")
  {
    double_param_ptr maxparam(new Parameter<double>(string(name+" Max"),pt.get_child("Max")));
    params["Max"] = maxparam;
    max_distance = params["Max"]->getLimHigh();
    distFunc = bind(&Uniform, _1, boost::ref(params["Max"]->getValue()));

  }
  else if (func_name == "Sine")
  {
    double_param_ptr maxparam(   new Parameter<double>(string(name+" Max"),    pt.get_child("Max")));
    double_param_ptr periodparam(new Parameter<double>(string(name+" Period"), pt.get_child("Period")));
    double_param_ptr offsetparam(new Parameter<double>(string(name+" Offset"), pt.get_child("Offset")));
    params["Max"]    = maxparam;
    params["Period"] = periodparam;
    params["Offset"] = offsetparam;
    max_distance = params["Max"]->getLimHigh();
    distFunc = bind(&Sine, _1, 
      boost::ref(params["Period"]->getValue()), 
      boost::ref(params["Offset"]->getValue()), 
      boost::ref(params["Max"]->getValue()));

  }
  else
  {
    stringstream err;
    err << "ERROR: read distance could not find function with name " << func_name << endl;
    error(err.str());
  }
  
  
}

void Distance::write(ptree& pt)
{
  ptree& distance_node = pt.add("Distance", "");
  
  distance_node.put("<xmlattr>.name",  name);
  distance_node.put("<xmlattr>.distfunc",  func_name);
  
  if      (func_name == "Linear")
  {
    ptree& max_node = distance_node.add("Max", "");
    params["Max"]->write(max_node, mode->getPrecision());
  } 
  else if (func_name == "Trapezoid")
  {
    ptree& A_node = distance_node.add("A", "");
    params["A"]->write(A_node, mode->getPrecision());
    ptree& B_node = distance_node.add("B", "");
    params["B"]->write(B_node, mode->getPrecision());
  }
  else if (func_name == "Uniform")
  {
    ptree& max_node = distance_node.add("Max", "");
    params["Max"]->write(max_node, mode->getPrecision());
  }
  else if (func_name == "Sine")
  {
    ptree& max_node    = distance_node.add("Max   ", "");
    ptree& period_node = distance_node.add("Period", "");
    ptree& offset_node = distance_node.add("Offset", "");
    
    params["Max"]->write(max_node, mode->getPrecision());
    params["Period"]->write(period_node, mode->getPrecision());
    params["Offset"]->write(offset_node, mode->getPrecision());
  }
}


double  Distance::getDistFunc(double distance)
{
  return distFunc(distance); 
}

double Distance::getMaxDistance()
{
  return max_distance;
}

void Distance::print(ostream& os)
{
  os << "Type: " << func_name << endl;
  
  typedef map<string, double_param_ptr> map_type;

  foreach_(map_type::value_type& myPair, params)
  {
    os << setw(6) << myPair.first << "  ";
    myPair.second->print(os);
  }
}



/************************    DistanceContainer    ************************************/


/*    Constructors    */

DistanceContainer::DistanceContainer() {}

DistanceContainer::DistanceContainer(ptree& pt) { add(pt); }


/*    Getters   */

void DistanceContainer::getParameters(param_ptr_vector& p)
{
  typedef map<string, distance_ptr>::iterator i_type;
  for (i_type i=distances.begin(); i != distances.end(); i++)
  {
    distance_ptr& distance = i->second;
    distance->getParameters(p);
  }
}

void DistanceContainer::getAllParameters(param_ptr_vector& p)
{
  typedef map<string, distance_ptr>::iterator i_type;
  for (i_type i=distances.begin(); i != distances.end(); i++)
  {
    distance_ptr& distance = i->second;
    distance->getAllParameters(p);
  }
}

/*    Setters   */


void DistanceContainer::add(ptree& pt)
{
  if (pt.count("Distances")) // we are at a parent of TFs
  {
    ptree& distances_node = pt.get_child("Distances");
    foreach_(ptree::value_type const& distance, distances_node)
    {
      if (distance.first == "Distance")
      {
        string name = distance.second.get<string>("<xmlattr>.name");
        distance_ptr dist(new Distance( (ptree&) distance.second));
        dist->setMode(mode);
        distances[name] = dist;
      }
    }
  }
  else if (pt.count("Distance")) // we are at the TFs node
  {
    foreach_(ptree::value_type const& distance, pt)
    {

      if (distance.first == "Distance")
      {
        string name = distance.second.get<string>("<xmlattr>.name");
        distance_ptr dist(new Distance( (ptree&) distance.second));
        dist->setMode(mode);
        distances[name] = dist;
      }
    }
  } else // this is a TF node
  {
    string name = pt.get<string>("<xmlattr>.name");
    distance_ptr dist(new Distance( (ptree&) pt));
    dist->setMode(mode);
    distances[name] = dist;  
  }
} 

void DistanceContainer::print(ostream& os)
{
  os << endl;
  typedef map<string, distance_ptr>::iterator i_type;
  for (i_type i=distances.begin(); i != distances.end(); i++)
  {
    os << i->first << endl;
    i->second->print(os);
  }
  os << endl;
}

void DistanceContainer::write(ptree& pt)
{
  ptree& distances_node = pt.add("Distances","");
  
  typedef map<string, distance_ptr>::iterator i_type;
  for (i_type i=distances.begin(); i != distances.end(); i++)
  {
    i->second->write(distances_node);
  }
}

distance_ptr DistanceContainer::getDistance(string name)
{
  if (distances.find(name) == distances.end())
  {
    stringstream err;
    err << "ERROR: could not find distance with name " << name << endl;
    error(err.str());
  }
  return distances[name];
}

double DistanceContainer::distFunc(string name, double distance)
{
  return distances[name]->getDistFunc(distance);
}
