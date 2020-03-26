/*********************************************************************************
*                                                                                *
*     coeffects.cpp                                                              *
*                                                                                *
*     Holds information about coactivation and corepression                      *
*                                                                                *
*********************************************************************************/

#include "coeffects.h"

#include <boost/foreach.hpp>
#define foreach_ BOOST_FOREACH

/**************************    Coeffect Class   *********************************/

/*    Constructors    */


Coeffect::Coeffect():
  efficiency(double_param_ptr(new Parameter<double>("CoeffectEff","CoeffectEff")))
{}

Coeffect::Coeffect(distances_ptr distances, ptree& pt):
  efficiency(double_param_ptr(new Parameter<double>("CoeffectEff","CoeffectEff")))
{ 
  this->distances = distances;
  read(pt); 
}


/*    Getters   */

void Coeffect::getParameters(param_ptr_vector& p)
{
  if (efficiency->isAnnealed())
    p.push_back(efficiency);
}

void Coeffect::getAllParameters(param_ptr_vector& p)
{
  p.push_back(efficiency);
}

  
/*    I/O   */

void Coeffect::read(ptree& pt)
{
  actor    = pt.get<string>("<xmlattr>.actor");
  target   = pt.get<string>("<xmlattr>.target");
  coef_idx = pt.get<int>("<xmlattr>.coef_idx");
  
  string distname = pt.get<string>("<xmlattr>.distance");
  dist = distances->getDistance(distname);
  
  ptree& orientation_node = pt.get_child("Orientation");
  
  HH = orientation_node.get<bool>("<xmlattr>.HH");
  HT = orientation_node.get<bool>("<xmlattr>.HT");
  TT = orientation_node.get<bool>("<xmlattr>.TT");
  
  efficiency->read(pt.get_child("efficiency"));
}

void Coeffect::write(ptree& pt)
{
  ptree& coeffect_node = pt.add("Coeffect","");
  
  coeffect_node.put("<xmlattr>.actor", actor);
  coeffect_node.put("<xmlattr>.target", target);
  coeffect_node.put("<xmlattr>.distance", dist->getName());
  coeffect_node.put("<xmlattr>.coef_idx", coef_idx);
  coeffect_node.put("<xmlattr>.include", true);
  
  ptree& efficiency_node = coeffect_node.add("efficiency", "");
  efficiency->write(efficiency_node, mode->getPrecision());
    
  ptree& orientation_node = coeffect_node.add("Orientation","");
  
  orientation_node.put("<xmlattr>.HH", HH);
  orientation_node.put("<xmlattr>.HT", HT);
  orientation_node.put("<xmlattr>.TT", TT);
}
  


/******************    Cooporativity Container Class   **************************/

/*    Constructors    */

CoeffectContainer::CoeffectContainer() {}

CoeffectContainer::CoeffectContainer(distances_ptr distances, ptree& pt) 
{
  read(pt, distances); 
}


/*    Getters   */


coeffect_pairs CoeffectContainer::getTargets(string tfname)
{
  coeffect_pairs output;
  int ncoeffects = coeffects.size();
  for (int i=0; i<ncoeffects; i++)
  {
    string& actor = coeffects[i]->getActor();

    if (actor == tfname)
    {
      string& target = coeffects[i]->getTarget();
    output.push_back(pair<string,coeffect_ptr>(target,coeffects[i]));
    }
  }
  return output;
}

void CoeffectContainer::getParameters(param_ptr_vector& p)
{
  int ncoefs = coeffects.size();
  for (int i=0; i<ncoefs; i++)
    coeffects[i]->getParameters(p);
}

void CoeffectContainer::getAllParameters(param_ptr_vector& p)
{
  int ncoefs = coeffects.size();
  for (int i=0; i<ncoefs; i++)
    coeffects[i]->getAllParameters(p);
}


/*    I/O   */

void CoeffectContainer::read(ptree& pt, distances_ptr distances)
{
  this->distances = distances;
  
  ptree& interactions_node = pt.get_child("Interactions");
  foreach_(ptree::value_type& node, interactions_node)
  {
    if (node.first != "Coeffect") continue;
    
    if ( !node.second.get<bool>("<xmlattr>.include", true) ) continue;
    
    coeffect_ptr cur_coef(new Coeffect(distances, (ptree&) node.second));
    cur_coef->setMode(mode);
    coeffects.push_back(cur_coef);
  }
}

void CoeffectContainer::write(ptree& pt)
{
  ptree* interactions_node;
  
  if (pt.count("Interactions") == 0)
    interactions_node = &(pt.add("Interactions", ""));
  else
    interactions_node = &(pt.get_child("Interactions"));
    
  int ncoefs = coeffects.size();
  for (int i=0; i<ncoefs; i++)
    coeffects[i]->write(*interactions_node);
}

/*
void CooperativityContainer::print(ostream& os)
{
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
    coops[i]->print(os);
}*/
