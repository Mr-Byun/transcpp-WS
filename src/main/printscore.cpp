/*********************************************************************************
*                                                                                *
*     unfold.cpp                                                                 *
*                                                                                *
*     Reads output and allows user to print various data                         *
*                                                                                *
*********************************************************************************/
#include "pwm.h"
#include "TF.h"
#include "gene.h"

#include "datatable.h"
#include "twobit.h"
#include "organism.h"

#include <fstream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>
#include <boost/lexical_cast.hpp>

#include <unistd.h>
#include <getopt.h>

using boost::property_tree::ptree;
#define to_string_ boost::lexical_cast<string>


int mode_verbose;
  
static const char *optString = "AcFfhMNorsXx:";

static const struct option longOpts[] = {
    { "help",         no_argument,       NULL, 'h' },
    { "section",      required_argument, NULL, 'x' },
    { "sequence",     no_argument,       NULL, 's' }
    //{ 0, 0, 0, 0}
};

void display_usage()
{
  cerr << endl << "\t Usage" << endl << endl
       << "\t printscore [options] input_file" << endl << endl
       << "\t Options" << endl
       << "\t --help       [-h]  print this message" << endl
       << "\t --section    [-x]  use section of input file (default eqparms)" << endl
       << "\t --sequence   [-s]  use sequence-level code" << endl << endl;
  exit(1);
}

int main(int argc, char* argv[])
{
  int opt = 0;
  int longIndex = 0;
  string section_name("eqparms");
  
  bool sequence = false;

  string infile_name;
  
  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while(opt != -1)
  {
    switch (opt)
    {
      case 'h':
        display_usage();
        break;
      case 's':
        sequence = true;
        break;
      case 'x':
        section_name = optarg;
        break;
      case 0: 
        display_usage();
        break;
      case '?':
        cerr << "?" << endl;
        display_usage();
        break;
    }
    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  }
  
  // convert to new node names
  if (section_name == "eqparms")
    section_name = "Output";
  else if (section_name == "input")
    section_name = "Input";
  
  if (argc <= optind)
    error("Missing input file");
  
  infile_name = argv[optind];
  
  if (infile_name == "")
    display_usage();
  
  ifstream infile(infile_name.c_str());

  ptree pt;
  read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);
  
  ptree& root_node   = pt.get_child("Root");
  ptree& mode_node   = root_node.get_child("Mode");
  ptree& section_node = root_node.get_child(section_name);
  
  mode_ptr mode(new Mode(infile_name,mode_node));
  
  mode->setScoreFunction("sse");
  mode->setScaleData(false);
  mode->setPerGene(false);
  mode->setPerNuc(false);

  mode->setVerbose(0);
  Organism embryo(section_node, mode);
  
  nuclei_ptr nuclei = embryo.getNuclei();
  genes_ptr  genes  = embryo.getGenes();
  
  int ngenes = genes->size();
  int nnuc   = nuclei->size();
  int ndata  = nnuc * ngenes;
  
  cout << " chisq = " << embryo.get_score() << endl;
  cout << " rms = " << sqrt(embryo.get_score()/ndata) << endl;
  return 0;
}


