/*********************************************************************************
*                                                                                *
*     scalefactor.cpp                                                            *
*                                                                                *
*     Contains a simple structure for scale factors                              *
*                                                                                *
*********************************************************************************/

#include "scalefactor.h"


#include <boost/foreach.hpp>
//#include <boost/range/adaptor/map.hpp>

#define foreach_ BOOST_FOREACH

/************************    ScaleFactor Class   ********************************/

/*    Constructors    */

ScaleFactor::ScaleFactor() 
{
  name = "default";
  double_param_ptr Aparam(new Parameter<double>("ScaleFactor","Null"));
  double_param_ptr Bparam(new Parameter<double>("ScaleFactor","Null"));
  A = Aparam;
  B = Bparam;
  A->set(1);
  B->set(0);
  A->setAnnealed(false);
  B->setAnnealed(false);
  A->setLimits(1,1);
  B->setLimits(0,0);
}

ScaleFactor::ScaleFactor(ptree& pt) 
{
  read(pt);
}


/*    Getters   */

void ScaleFactor::getParameters(param_ptr_vector& p)
{
  if (A->isAnnealed())
    p.push_back(A);
  if (B->isAnnealed())
    p.push_back(B);
}

void ScaleFactor::getAllParameters(param_ptr_vector& p)
{
  p.push_back(A);
  p.push_back(B);
}


/*    Methods   */

double ScaleFactor::scale(double x)
{
  return max(x * A->getValue() + B->getValue(), 0.0);
}

double ScaleFactor::unscale(double x)
{
  return max( (x - B->getValue())/A->getValue(), 0.0);
}


/*    I/O   */

void ScaleFactor::read(ptree& pt)
{
  name = pt.get<string>("<xmlattr>.name");
  double_param_ptr Aparam(new Parameter<double>(pt.get_child("A"),"ScaleFactor","Null"));
  double_param_ptr Bparam(new Parameter<double>(pt.get_child("B"),"ScaleFactor","Null"));
  A = Aparam;
  B = Bparam;
}

void ScaleFactor::write(ptree& pt)
{
  ptree& factor_node = pt.add("ScaleFactor", "");
  factor_node.put("<xmlattr>.name", name);
  ptree& A_node = factor_node.add("A", "");
  A->write(A_node, mode->getPrecision());
  ptree& B_node = factor_node.add("B", "");
  B->write(B_node, mode->getPrecision());
}


/*******************    ScaleFactor Container Class   ***************************/


/*    Constructors    */

ScaleFactorContainer::ScaleFactorContainer() {}

ScaleFactorContainer::ScaleFactorContainer(ptree& pt) 
{
  read(pt);
}


/*    Getters   */

void ScaleFactorContainer::getParameters(param_ptr_vector& p)
{
  //foreach_(scale_factor_ptr scale, scales | map_values)
  //  scale->getParameters(p);
  int nscales = scales_vector.size();
  for (int i=0; i<nscales; i++)
    scales_vector[i]->getParameters(p);
}

void ScaleFactorContainer::getAllParameters(param_ptr_vector& p)
{
  //foreach_(scale_factor_ptr scale, scales | map_values)
  //  scale->getAllParameters(p);
  int nscales = scales_vector.size();
  for (int i=0; i<nscales; i++)
    scales_vector[i]->getAllParameters(p);
}

scale_factor_ptr ScaleFactorContainer::getScaleFactor(string name) 
{ 
  if (scales.find(name) == scales.end())
    return scales["default"];
  else
    return scales[name];
}

/*    I/O   */

void ScaleFactorContainer::read(ptree& pt)
{
  ptree& scale_factors_node = pt.get_child("ScaleFactors");
  
  foreach_(ptree::value_type& node, scale_factors_node)
  {
    if (node.first != "ScaleFactor") continue;
    
    scale_factor_ptr scale(new ScaleFactor((ptree&) node.second));
    scale->setMode(mode);
    string& name = scale->getName();
    scales[name] = scale;
    scales_vector.push_back(scale);
  }
  
  scale_factor_ptr scale(new ScaleFactor());
  scale->setMode(mode);
  scales["default"] = scale;
  scales_vector.push_back(scale);
}

void ScaleFactorContainer::write(ptree& pt)
{
  ptree& scales_node = pt.add("ScaleFactors", "");
  
  //foreach_(scale_factor_ptr scale, scales | map_values)
  //  scale->write(scales_node);

  int nscales = scales_vector.size();
  for (int i=0; i<nscales; i++)
    scales_vector[i]->write(scales_node);
}
