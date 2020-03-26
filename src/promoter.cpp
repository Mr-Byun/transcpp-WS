/*********************************************************************************
*                                                                                *
*     promoter.cpp                                                               *
*                                                                                *
*     Contains parameters that are shared among many nuclei (Q,theta)            *
*                                                                                *
*********************************************************************************/

#include "promoter.h"
#include "utils.h"

#include <math.h>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <boost/ref.hpp>
#include <limits>

# define foreach_ BOOST_FOREACH


/************************    Functions   ****************************************/

// the function in the old code. Less stable for large M values
double Arrhenius(double M, double max, double theta, double Q)
{
  double exponent = exp(Q*M - theta);
  return max * (exponent) / (1 + exponent);
}

// Same eq as 1, but simpler and more stable
// same output as arrhenius, but never returns NaN
double Arrhenius2(double M, double max, double theta, double Q)
{
  double exponent = exp(theta - Q*M);
  return max / (1 + exponent);
}

// the old exponential rate
double Exponential(double M, double max, double theta, double Q)
{
  double exponent = exp(Q*M - theta);
  if (exponent > max)
    return max;
  else
    return exponent;
}

double Linear(double M, double A, double B)
{
  return M*A + B;
}


/************************    Promoter Class   ***********************************/

/*    Constructors    */

Promoter::Promoter() {}

Promoter::Promoter(ptree& pt) { read(pt); }


/*    Getters   */

void Promoter::getParameters(param_ptr_vector& p)
{
  typedef map<string, double_param_ptr>::iterator i_type;
  for(i_type i=params.begin(); i != params.end(); i++)
  {
    double_param_ptr& param = i->second;
    if (param->isAnnealed())
      p.push_back(param);
  }
}

void Promoter::getAllParameters(param_ptr_vector& p)
{
  typedef map<string, double_param_ptr>::iterator i_type;
  for(i_type i=params.begin(); i != params.end(); i++)
  {
    double_param_ptr& param = i->second;
    p.push_back(param);
  }
}


/*    I/O   */

void Promoter::read(ptree& pt)
{
  name      = pt.get<string>("<xmlattr>.name");
  func_name = pt.get<string>("<xmlattr>.function");
  
  if (func_name == "Arrhenius")
  {
    double_param_ptr qparam(    new Parameter<double>(pt.get_child("Q"), string(name+" Q"),"Promoter"));
    double_param_ptr maxparam(  new Parameter<double>(pt.get_child("Rmax"),string(name+" Rmax"),"Promoter"));
    double_param_ptr thetaparam(new Parameter<double>(pt.get_child("Theta"),string(name+" Theta"),"Promoter"));
    
    params["Q"]     = qparam;
    params["Rmax"]  = maxparam;
    params["Theta"] = thetaparam;
    
    rateFunc = bind(&Arrhenius, _1, 
      boost::ref(params["Rmax"]->getValue()), 
      boost::ref(params["Theta"]->getValue()), 
      boost::ref(params["Q"]->getValue()));
  }
  else if (func_name == "Arrhenius2")
  {
    double_param_ptr qparam(    new Parameter<double>(pt.get_child("Q"), string(name+" Q"),"Promoter"));
    double_param_ptr maxparam(  new Parameter<double>(pt.get_child("Rmax"),string(name+" Rmax"),"Promoter"));
    double_param_ptr thetaparam(new Parameter<double>(pt.get_child("Theta"),string(name+" Theta"),"Promoter"));
    
    params["Q"]     = qparam;
    params["Rmax"]  = maxparam;
    params["Theta"] = thetaparam;
    
    rateFunc = bind(&Arrhenius2, _1, 
      boost::ref(params["Rmax"]->getValue()), 
      boost::ref(params["Theta"]->getValue()), 
      boost::ref(params["Q"]->getValue()));
  }
  else if (func_name == "Exponential")
  {
    double_param_ptr qparam(    new Parameter<double>(pt.get_child("Q"), string(name+" Q"),"Promoter"));
    double_param_ptr maxparam(  new Parameter<double>(pt.get_child("Rmax"),string(name+" Rmax"),"Promoter"));
    double_param_ptr thetaparam(new Parameter<double>(pt.get_child("Theta"),string(name+" Theta"),"Promoter"));
    
    params["Q"]     = qparam;
    params["Rmax"]  = maxparam;
    params["Theta"] = thetaparam;
    
    rateFunc = bind(&Exponential, _1, 
      boost::ref(params["Rmax"]->getValue()), 
      boost::ref(params["Theta"]->getValue()), 
      boost::ref(params["Q"]->getValue()));
    
  }
  else if (func_name == "Linear")
  {
    double_param_ptr Aparam(  new Parameter<double>(pt.get_child("A"),"A","Promoter"));
    double_param_ptr Bparam(  new Parameter<double>(pt.get_child("B"),"A","Promoter"));
    
    params["A"]  = Aparam;
    params["B"]  = Bparam;
    
    rateFunc = bind(&Linear, _1, 
      boost::ref(params["A"]->getValue()), 
      boost::ref(params["B"]->getValue()));
  }
  else
  {
    stringstream err;
    err << "ERROR: read promoter could not find function with name " << func_name << endl;
    error(err.str());
  }
}
    
void Promoter::write(ptree& pt)
{
  ptree& promoter_node = pt.add("Promoter","");
  
  promoter_node.put("<xmlattr>.name",     name);
  promoter_node.put("<xmlattr>.function", func_name);
  
  if (func_name == "Arrhenius")
  {
    ptree& q_node     = promoter_node.add("Q    ","");
    ptree& max_node   = promoter_node.add("Rmax ","");
    ptree& theta_node = promoter_node.add("Theta","");
    
    params["Q"]->write(q_node, mode->getPrecision());
    params["Rmax"]->write(max_node, mode->getPrecision());
    params["Theta"]->write(theta_node, mode->getPrecision());
  }
  else if (func_name == "Arrhenius2")
  {
    ptree& q_node     = promoter_node.add("Q    ","");
    ptree& max_node   = promoter_node.add("Rmax ","");
    ptree& theta_node = promoter_node.add("Theta","");
    
    params["Q"]->write(q_node, mode->getPrecision());
    params["Rmax"]->write(max_node, mode->getPrecision());
    params["Theta"]->write(theta_node, mode->getPrecision());
  }
  else if (func_name == "Exponential")
  {
    ptree& q_node     = promoter_node.add("Q    ","");
    ptree& max_node   = promoter_node.add("Rmax ","");
    ptree& theta_node = promoter_node.add("Theta","");
    
    params["Q"]->write(q_node, mode->getPrecision());
    params["Rmax"]->write(max_node, mode->getPrecision());
    params["Theta"]->write(theta_node, mode->getPrecision());
  }
  else if (func_name == "Linear")
  {
    ptree& A_node   = promoter_node.add("A ","");
    ptree& B_node   = promoter_node.add("B ","");
    
    params["A"]->write(A_node, mode->getPrecision());
    params["B"]->write(B_node, mode->getPrecision());

  } else
    error("Could not find write function for " + func_name);
}


/*******************    Promoter Container Class   ******************************/


/*    Constructors    */

PromoterContainer::PromoterContainer() {}

PromoterContainer::PromoterContainer(ptree& pt) { read(pt); }


/*    Getters   */

void PromoterContainer::getParameters(param_ptr_vector& p)
{
  typedef map<string, promoter_ptr>::iterator i_type;
  for(i_type i=promoters.begin(); i != promoters.end(); i++)
  {
    promoter_ptr& promoter = i->second;
    promoter->getParameters(p);
  }
}

void PromoterContainer::getAllParameters(param_ptr_vector& p)
{
  typedef map<string, promoter_ptr>::iterator i_type;
  for(i_type i=promoters.begin(); i != promoters.end(); i++)
  {
    promoter_ptr& promoter = i->second;
    promoter->getAllParameters(p);
  }
}


/*    I/O   */

void PromoterContainer::read(ptree& pt)
{
  ptree& promoters_node = pt.get_child("Promoters");
  
  foreach_(ptree::value_type& node, promoters_node)
  {
    if (node.first != "Promoter") continue;
    
    promoter_ptr promoter(new Promoter((ptree&) node.second));
    promoter->setMode(mode);
    string& name = promoter->getName();
    promoters[name] = promoter;
  }
}

void PromoterContainer::write(ptree& pt)
{
  ptree& promoters_node = pt.add("Promoters", "");
  
  typedef map<string, promoter_ptr>::iterator i_type;
  for(i_type i=promoters.begin(); i != promoters.end(); i++)
  {
    promoter_ptr& promoter = i->second;
    promoter->write(promoters_node);
  }
}
    
