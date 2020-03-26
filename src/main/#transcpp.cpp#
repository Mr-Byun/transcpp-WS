/*********************************************************************************
*                                                                                *
*     transc.cpp                                                                 *
*                                                                                *
*     Constains the main method for fitting the transcription model              *
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
#include <boost/version.hpp>

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


  string xmlname(argv[1]);
  fstream infile(xmlname.c_str());
  
  ptree pt;
  read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);
  
  ptree& root_node  = pt.get_child("Root");
  ptree& mode_node  = root_node.get_child("Mode");
  ptree& input_node = root_node.get_child("Input");
  
  mode_ptr mode(new Mode(xmlname, mode_node));
  
  if (mode->getVerbose() >= 1)
  {
    cerr << endl;
    cerr << "Input file read" << endl;
    cerr << "Verbose level set to "      << mode->getVerbose() << endl;
    cerr << endl;
    if (mode->getSchedule() == LAM)
      cerr << "Annealing schedule set to " << mode->getSchedule() << ": lam." << endl;
    else if (mode->getSchedule() == EXP)
      cerr << "Annealing schedule set to " << mode->getSchedule() << ": exp." << endl;
    else
      error("Could not identify annealing schedule");
    cerr << endl;
  }

  Organism embryo(input_node, mode);
  //embryo.printParameters(cerr);
  
  unirand48 rnd;
  unsigned int seed = mode->getSeed();
  if (mode->getVerbose() >= 1) cerr << "Beginning annealing with seed " << seed << endl;
  rnd.setSeed(seed);
  
  xmlDoc *doc = xmlParseFile(xmlname.c_str());
  xmlNode *docroot = xmlDocGetRootElement(doc);

  criCount::Param criCntParam(docroot);
  
  annealer<Organism, lam, criCount, feedbackMove>*     fly_sa;
  annealer<Organism, expHold, criCount, feedbackMove>* fly_expHold;
  
  if (mode->getSchedule() == LAM)
  {
    lam::Param scheduleParam(docroot);
    fly_sa = new annealer<Organism, lam, criCount, feedbackMove>(embryo, rnd, scheduleParam, criCntParam, docroot);
    fly_sa->setCoolLog(file, (xmlname+".log").c_str());
    fly_sa->setProlix(file, (xmlname+".prolix").c_str());
    if (mode->getVerbose() >= 1)
      cerr << "The initial score is " << embryo.get_score() << endl;
    fly_sa->initMoves();
    if (mode->getVerbose() >= 1)
      cerr << "The score is " << embryo.get_score() << " after initial moves" << endl << endl;
    if (!mode->getProfiling())
    {
      fly_sa->loop();
      if (mode->getVerbose() >= 1)
      {
        cerr << "The final score is " << embryo.get_score() << endl << endl;
        fly_sa->writeResult();
      }
    }
  } 
  else if (mode->getSchedule() == EXP)
  {
    expHold::Param scheduleParam(docroot);
    fly_expHold = new annealer<Organism, expHold, criCount, feedbackMove>(embryo, rnd, scheduleParam, criCntParam, docroot);
    fly_expHold->setStepLog(file, (xmlname+".steplog").c_str());
    fly_expHold->setProlix(file,(xmlname+".prolix").c_str());
    if (!mode->getProfiling())
      fly_expHold->loop();
  } 


  	
  if (!mode->getProfiling())
  {
    embryo.write("Output", root_node);

#if BOOST_VERSION / 100 % 1000 < 56
    write_xml(xmlname, 
      pt, 
      std::locale(), 
      boost::property_tree::xml_writer_make_settings<char>(' ', 2));
#else
    write_xml(xmlname, 
      pt, 
      std::locale(), 
      boost::property_tree::xml_writer_make_settings<string>(' ', 2));
#endif
    

    /* check that I can reset everything and get the same correct answer */
    double old_score = embryo.get_score();
    embryo.ResetAll();
    double score     = embryo.get_score();
    
    if (score != old_score)
      error("The score was not the same after ResetAll. Something is not being moved properly!");
    
    /* here I double check the output and make sure it gives me the right score*/
    
    if (mode->getVerbose() >= 2)
      cerr << endl << "verifying output..." << endl;
    score = embryo.get_score();
    
    infile.close();
    
    fstream outfile(xmlname.c_str());
  
    ptree pt_out;
    read_xml(outfile, pt_out, boost::property_tree::xml_parser::trim_whitespace);
  
    ptree& root_node_out  = pt_out.get_child("Root");
    ptree& mode_node_out  = root_node_out.get_child("Mode");
    ptree& input_node_out = root_node_out.get_child("Output");
    
    mode_ptr mode_out(new Mode(xmlname, mode_node_out));

    Organism embryo_out(input_node_out, mode_out);

    //embryo_out.printRate(cerr, 0);
    
    /* this has to be done with some precision because we arent always using 16 bit.
    I am going to round to half the decimal places in precision */
    double digits = mode->getPrecision()/2 * 10;
    double score_out = embryo_out.get_score();
    
    score     = round(score*digits)/digits;
    score_out = round(score_out*digits)/digits;
    
    if (score != score_out)
      error("The score (" + to_string_(score) + 
        ") was not the same after rerunning the model (" +
        to_string_(score_out) +"). Something may not have been printed!");
    else
    {
      if (mode->getVerbose()>=2)
        cerr << "output verified!" << endl << endl;
    }
    
  }
  
  xmlFreeDoc(doc);
  
  return 0;
}
