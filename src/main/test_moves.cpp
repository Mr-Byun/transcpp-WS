/*********************************************************************************
*                                                                                *
*     test_moves.cpp                                                             *
*                                                                                *
*     This function tests whether or not annealing appears to be working         *
*     for a given input file my moving everything individually and making sure   *
*     that reset all and recalculate all give the same answer                    *
*                                                                                *
*********************************************************************************/

#include "pwm.h"
#include "TF.h"
#include "gene.h"
#include "datatable.h"
#include "twobit.h"
#include "organism.h"
#include "mode.h"
#include "utils.h"
#include <fstream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <unistd.h>
#include <libxml/parser.h>
#include "annealer.h"
#include "move/feedbackMove.h"
#include "unirandom.h"
#include "lam.h"
#include "expHold.h"
#include "tempCount.h"
#include "criCount.h"
#include "dynDebug.h"

#define to_        boost::lexical_cast
#define to_string_ boost::lexical_cast<string>


using boost::property_tree::ptree;

int main(int argc, char* argv[])
{
  /* create my random number generator */
  boost::minstd_rand baseGen(getpid());
  boost::uniform_real<> uniDblUnit(0,1);
  boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uniDblGen(baseGen, uniDblUnit);
  uniDblGen();
  

  string xmlname(argv[1]);
  fstream infile(xmlname.c_str());
  
  ptree pt;
  read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);
  
  ptree& root_node  = pt.get_child("Root");
  ptree& mode_node  = root_node.get_child("Mode");
  ptree& input_node = root_node.get_child("Input");
  
  mode_ptr mode(new Mode(xmlname, mode_node));
  
  mode->setVerbose(0);
  for (int i=0; i<10; i++)
    Organism embryo(input_node, mode);
  
  Organism embryo(input_node, mode);
  
  unirand48 rnd;
  //unsigned int seed = mode->getSeed();
  
  cerr << endl;
  
  int nparams = embryo.getDimension();
  for (int i=0; i<nparams; i++)
  {
    iparam_ptr iparam = embryo.getParam(i);
    ParameterInterface* param = iparam.get();
    cerr << "Testing parameter " << i << ": " << embryo.getParamName(i) << endl;
    
    if (param->getType() == "PWM")
    {
      for (int j=0; j<10; j++)
      {
        double rand       = uniDblGen();
        double delta    = 10*rand - 5;
        double start_score = embryo.get_score();
        
        embryo.generateMove(i, delta);
        while(embryo.get_score() == numeric_limits<double>::max())
          embryo.generateMove(i, delta);
        
        
        embryo.restoreMove(i);
        
       
        if (embryo.get_score() != start_score)
          error("Starting score was " + to_string_(start_score) + " but restore gave " + to_string_(embryo.get_score()));

        embryo.generateMove(i, delta);
        while(embryo.get_score() == numeric_limits<double>::max())
          embryo.generateMove(i, delta);
        double  move_score      = embryo.get_score();
        //cerr << "move: " << move_score << endl;
        //Bindings move_bindings = *(embryo.getBindings());
        embryo.ResetAll();
        double reset_score  = embryo.get_score();
        //cerr << "reset: " << reset_score << endl;
        //Bindings reset_bindings = *(embryo.getBindings());
        embryo.Recalculate();
        double recalc_score = embryo.get_score();
        //cerr << "recalc: " << recalc_score << endl;
        //Bindings recalc_bindings = *(embryo.getBindings());
        
        if (move_score != reset_score)
        {
          //move_bindings.isEqual(reset_bindings);
          error("Move ("+ to_string_(move_score)+") and ResetAll ("+ to_string_(reset_score)+") gave different answers. Move generation may be broken");
        }
        if (reset_score != recalc_score)
        {
          //reset_bindings.isEqual(recalc_bindings);
          error("ResetAll ("+ to_string_(reset_score)+") and Recalculate ("+ to_string_(recalc_score)+") gave different answers. ResetAll may be broken");
        }
        cerr << ".";
      }
      cerr << endl;
    }
    else
    {
      
      Parameter<double>* p = dynamic_cast<Parameter<double>* >(param);
      
      // generate a delta that will be in bounds
      double lim_low  = p->getLimLow();
      double lim_high = p->getLimHigh();
      
      
      for (int j=0; j<10; j++)
      {
        if (param->getType() != string("double")) continue;
        double value    = p->getValue();
        double pmax     = lim_high - value;
        double pmin     = lim_low - value;
        
        //cerr << "range: " << pmin << "-" << pmax << endl;
        double start_score = embryo.get_score();
        
        double rand     = uniDblGen();
        double delta    = (pmax - pmin)*rand + pmin;
        
        //cerr << "value: " << value << endl;
        //cerr << "delta: " << delta << endl;
        
        
        // check to see if restore is working
        embryo.generateMove(i, delta);
        if (p->isOutOfBounds())
          error("Parameter ("+to_string_(value+delta)+") is out of bounds ("+to_string_(lim_low)+","+to_string_(lim_high)+"). Check move generation in this file!");
        embryo.restoreMove(i);
        if (embryo.get_score() != start_score)
          error("Restore Failed!");
        
        // check to see that the move gives the same score as reset all and recalculate
        embryo.generateMove(i, delta);
        if (p->isOutOfBounds())
          error("Parameter ("+to_string_(value+delta)+") is out of bounds ("+to_string_(lim_low)+","+to_string_(lim_high)+"). Check move generation in this file!");
        double move_score      = embryo.get_score();
        //cerr << "move: " << move_score << endl;
        //Bindings move_bindings = *(embryo.getBindings());
        embryo.ResetAll();
        double reset_score  = embryo.get_score();
        //cerr << "reset: " << reset_score << endl;
        //Bindings reset_bindings = *(embryo.getBindings());
        embryo.Recalculate();
        double recalc_score = embryo.get_score();
        //cerr << "recalc: " << recalc_score << endl;
        //Bindings recalc_bindings = *(embryo.getBindings());
       
        if (move_score != reset_score)
        {
          //move_bindings.isEqual(reset_bindings);
          error("Move ("+ to_string_(move_score)+") and ResetAll ("+ to_string_(reset_score)+") gave different answers. Move generation may be broken");
        }
        if (reset_score != recalc_score)
        {
          //reset_bindings.isEqual(recalc_bindings);
          error("ResetAll ("+ to_string_(reset_score)+") and Recalculate ("+ to_string_(recalc_score)+") gave different answers. ResetAll may be broken");
        }
        
        cerr << ".";
      }
      cerr << endl;
    }
  }
  
  cerr << endl << "All move functions appear to be working for this problem! Congratulations!" << endl << endl;
}
  
  
