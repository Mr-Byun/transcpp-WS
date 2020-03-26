/*********************************************************************************
*                                                                                *
*     competition.cpp                                                            *
*                                                                                *
*     A very simple class holding the parameters for promoter competition        *
*                                                                                *
*********************************************************************************/

#include "competition.h"
#include "utils.h"
#include "mode.h"


/*    Constructors    */

Competition::Competition():
  window(double_param_ptr(new Parameter<double>("Window","Window"))),
  shift(double_param_ptr(new Parameter<double>("Shift","Promoter"))),
  specificity(double_param_ptr(new Parameter<double>("Specificity","Promoter"))),
  threshold(double_param_ptr(new Parameter<double>("CompThreshold","Promoter"))),
  background(double_param_ptr(new Parameter<double>("CompBackground","Promoter"))),
  S(double_param_ptr(new Parameter<double>("S","Promoter"))),
  interaction_strength(double_param_ptr(new Parameter<double>("InteractionStrength","Promoter")))
{}

Competition::Competition(ptree& pt, mode_ptr mode):
  window(double_param_ptr(new Parameter<double>("Window","Window"))),
  shift(double_param_ptr(new Parameter<double>("Shift","Promoter"))),
  specificity(double_param_ptr(new Parameter<double>("Specificity","Promoter"))),
  threshold(double_param_ptr(new Parameter<double>("CompThreshold","Promoter"))),
  background(double_param_ptr(new Parameter<double>("CompBackground","Promoter"))),
  S(double_param_ptr(new Parameter<double>("S","Promoter"))),
  interaction_strength(double_param_ptr(new Parameter<double>("InteractionStrength","Promoter")))
{
  set(pt, mode); 
}


/*    Getters   */

void Competition::getParameters(param_ptr_vector& p)
{
  if (window->isAnnealed())
    p.push_back(window);
  if (shift->isAnnealed())
    p.push_back(shift);
  if (specificity->isAnnealed())
    p.push_back(specificity);
  if (threshold->isAnnealed())
    p.push_back(threshold);
  if (background->isAnnealed())
    p.push_back(background);
  if (S->isAnnealed())
    p.push_back(S);
  if (interaction_strength->isAnnealed())
    p.push_back(interaction_strength);
}

void Competition::getAllParameters(param_ptr_vector& p)
{
  p.push_back(window);
  p.push_back(shift);
  p.push_back(specificity);
  p.push_back(threshold);
  p.push_back(background);
  p.push_back(S);
  p.push_back(interaction_strength);
}


/*    Setters   */


/*    I/O   */

void Competition::set(ptree& pt, mode_ptr mode)
{
  this->mode = mode;
  if (pt.count("Competition") == 0)
  {
    warning("no competition node: setting values in mode");
    window->set(mode->getWindow());
    shift->set(mode->getShift());
    specificity->set(mode->getN());
    threshold->set(mode->getT());
    interaction_strength->set(1);
    background->set(0);
    S->set(1);
  }
  else
  {
    ptree& competition_node = pt.get_child("Competition");
    
    if (competition_node.count("Window") == 0)
      window->set(mode->getWindow());
    else
    {
      ptree& window_node = competition_node.get_child("Window");
      window->read(window_node);
    }
    
    if (competition_node.count("Shift") == 0)
      shift->set(mode->getShift());
    else
    {
      ptree& shift_node = competition_node.get_child("Shift");
      shift->read(shift_node);
    }
    
    if (competition_node.count("Specificity") == 0)
      specificity->set(mode->getN());
    else
    {
      ptree& specificity_node = competition_node.get_child("Specificity");
      specificity->read(specificity_node);
    }
    
    if (competition_node.count("Threshold") == 0)
      threshold->set(mode->getT());
    else
    {
      ptree& threshold_node = competition_node.get_child("Threshold");
      threshold->read(threshold_node);
    }
    
    if (competition_node.count("Background") == 0)
      background->set(0);
    else
    {
      ptree& background_node = competition_node.get_child("Background");
      background->read(background_node);
    }
    
    if (competition_node.count("S") == 0)
      S->set(1);
    else
    {
      ptree& S_node = competition_node.get_child("S");
      S->read(S_node);
    }
    
    if (competition_node.count("InteractionStrength") == 0)
      interaction_strength->set(1);
    else
    {
      ptree& strength_node = competition_node.get_child("InteractionStrength");
      interaction_strength->read(strength_node);
    }
    
    if (competition_node.count("NProportionality") == 0)
      product = false;
    else
    {
      ptree& nprop_node = competition_node.get_child("NProportionality");
      string str = nprop_node.get<string>("<xmlattr>.value", "sum");
      if (str == string("sum"))
        product = false;
      else if (str == string("product"))
        product = true;
      else
        error("Unrecognized proportionality type " + str);
    }
  }
}


void Competition::write(ptree& pt)
{
  ptree& competition_node = pt.add("Competition", "");
  
  ptree& window_node      = competition_node.add("Window     ", "");
  window->write(window_node, mode->getPrecision());
  
  ptree& shift_node       = competition_node.add("Shift      ", "");
  shift->write(shift_node, mode->getPrecision());
  
  ptree& specificity_node = competition_node.add("Specificity", "");
  specificity->write(specificity_node, mode->getPrecision());
  
  ptree& threshold_node   = competition_node.add("Threshold  ", "");
  threshold->write(threshold_node, mode->getPrecision());

  ptree& background_node  = competition_node.add("Background ", "");
  background->write(background_node, mode->getPrecision());
  
  ptree& S_node           = competition_node.add("S          ", "");
  S->write(S_node, mode->getPrecision());
  
  ptree& strength_node    = competition_node.add("InteractionStrength","");
  interaction_strength->write(strength_node, mode->getPrecision());
  
  ptree& nprop_node  = competition_node.add("NProportionality ", "");
  if (product)
    nprop_node.put("<xmlattr>.value", "product");
  else
    nprop_node.put("<xmlattr>.value", "sum");
}


  
  
