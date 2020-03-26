/*********************************************************************************
*                                                                                *
*     cooperativity.cpp                                                          *
*                                                                                *
*     Holds information about cooporativity between tfs                          *
*                                                                                *
*********************************************************************************/

#include "cooperativity.h"

#include <boost/foreach.hpp>
#define foreach_ BOOST_FOREACH

/***********************    Cooporativity Class   *******************************/

/*    Constructors    */


Cooperativity::Cooperativity():
  Kcoop(double_param_ptr(new Parameter<double>("Kcoop","Kcoop")))
{}

Cooperativity::Cooperativity(distances_ptr distances, ptree& pt):
  Kcoop(double_param_ptr(new Parameter<double>("Kcoop","Kcoop")))
{ 
  this->distances = distances;
  read(pt); 
}


/*    Getters   */

void Cooperativity::getParameters(param_ptr_vector& p)
{
  if (Kcoop->isAnnealed())
    p.push_back(Kcoop);
}

void Cooperativity::getAllParameters(param_ptr_vector& p)
{
  p.push_back(Kcoop);
}

pair<string,string> Cooperativity::getTFs()
{
  return pair<string, string>(factor1,factor2);
}
  
/*    I/O   */

void Cooperativity::read(ptree& pt)
{
  factor1 = pt.get<string>("<xmlattr>.factor1");
  factor2 = pt.get<string>("<xmlattr>.factor2");
  
  ptree& orientation_node = pt.get_child("Orientation");
  
  HH = orientation_node.get<bool>("<xmlattr>.HH");
  HT = orientation_node.get<bool>("<xmlattr>.HT");
  TT = orientation_node.get<bool>("<xmlattr>.TT");
  
  Kcoop->read(pt.get_child("Kcoop"));
  
  string distname = pt.get<string>("<xmlattr>.distance");
  dist = distances->getDistance(distname);
}

void Cooperativity::write(ptree& pt)
{
  ptree& coop_node = pt.add("Cooperativity","");
  coop_node.put("<xmlattr>.factor1", factor1);
  coop_node.put("<xmlattr>.factor2", factor2);
  
  ptree& orientation_node = coop_node.add("Orientation","");
  
  orientation_node.put("<xmlattr>.HH", HH);
  orientation_node.put("<xmlattr>.HT", HT);
  orientation_node.put("<xmlattr>.TT", TT);
  
  ptree& kcoop_node = coop_node.add("Kcoop", "");
  Kcoop->write(kcoop_node, mode->getPrecision());
  
  coop_node.put("<xmlattr>.distance", dist->getName());
}

void Cooperativity::print(ostream& os)
{
  os << endl
     << "\t Cooperativity " << endl
     << "\tfactor 1: " << factor1 << endl
     << "\tfactor 2: " << factor2 << endl
     << "\tHH: " << HH << "  HT: " << HT << "  TT: " << TT << endl;
  Kcoop->print(os);
  os << endl;
}

/******************    Cooporativity Container Class   **************************/

/*    Constructors    */

CooperativityContainer::CooperativityContainer() {}

CooperativityContainer::CooperativityContainer(distances_ptr distances, ptree& pt) 
{
  read(pt, distances); 
}


/*    Getters   */

coop_pairs CooperativityContainer::getCoops(string tfname)
{
  coop_pairs output;
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
  {
    pair<string,string> tf_pair = coops[i]->getTFs();
    if (tf_pair.first == tfname)
      output.push_back(pair<string,coop_ptr>(tf_pair.second,coops[i]));
    else if (tf_pair.second == tfname)
      output.push_back(pair<string,coop_ptr>(tf_pair.first,coops[i]));
  }
  return output;
}

void CooperativityContainer::getParameters(param_ptr_vector& p)
{
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
    coops[i]->getParameters(p);
}

void CooperativityContainer::getAllParameters(param_ptr_vector& p)
{
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
    coops[i]->getAllParameters(p);
}


/*    I/O   */

void CooperativityContainer::read(ptree& pt, distances_ptr distances)
{
  this->distances = distances;
  
  ptree& interactions_node = pt.get_child("Interactions");
  foreach_(ptree::value_type& node, interactions_node)
  {
    if (node.first != "Cooperativity") continue;
    
    if ( !node.second.get<bool>("<xmlattr>.include", true) ) continue;
    
    coop_ptr cur_coop(new Cooperativity(distances, (ptree&) node.second));
    cur_coop->setMode(mode);
    coops.push_back(cur_coop);
  }
}

void CooperativityContainer::write(ptree& pt)
{
  ptree* interactions_node;
  
  if (pt.count("Interactions") == 0)
    interactions_node = &(pt.add("Interactions", ""));
  else
    interactions_node = &(pt.get_child("Interactions"));
  
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
    coops[i]->write(*interactions_node);
}

void CooperativityContainer::print(ostream& os)
{
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
    coops[i]->print(os);
}
