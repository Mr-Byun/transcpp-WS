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

#include <unistd.h>
#include <getopt.h>

using boost::property_tree::ptree;


int mode_verbose;

static const char *optString = "hi:s:";

static const struct option longOpts[] = {
    { "help",        no_argument,       NULL, 'h' },
    { "input-file",  required_argument, NULL, 'i' },
    { "section",     required_argument, NULL, 's' },
    { "sites",       no_argument,       NULL,  0  },
    { "occupancy",   no_argument,       NULL,  0  },
    { "modeocc",     no_argument,       NULL,  0  },
    { "effocc",      no_argument,       NULL,  0  },
    { "subgroups",   no_argument,       NULL,  0  },
    { "scores",      no_argument,       NULL,  0  },
    { "section",     required_argument, NULL, 's' },
    { "gene",        required_argument, NULL,  0  },
    { "tf",          required_argument, NULL,  0  },
    { "rate",        no_argument,       NULL,  0  },
    { "R2D",         no_argument,       NULL,  0  },
    { "N2D",         no_argument,       NULL,  0  },
    { "data",        no_argument,       NULL,  0  },
    { "params",      no_argument,       NULL,  0  },
    { "check-scale", no_argument,       NULL,  0  },
    { "invert",      no_argument,       NULL,  0  },
    { "score",       no_argument,       NULL,  0  },
    { "maxscore",    no_argument,       NULL,  0  },
    { 0, 0, 0, 0}
};

void display_usage()
{
  cerr << endl << "\t Usage" << endl << endl
       << "\t unfold [options] -i [infile]" << endl << endl
       << "\t Options" << endl
       << "\t --help    [-h]   print this message" << endl
       << "\t --section [-s]   use section of input file (default Output)" << endl
       << "\t --score          prints the score of the entire fit" << endl
       << "\t --maxscore       prints max score of each TF" << endl
       << "\t --sites          prints binding sites" << endl
       << "\t --occupancy      prints binding site fractional occupancy" << endl
       << "\t --modeocc        prints the occupancy split into tf modes" << endl
       << "\t --effocc         prints the effective (activating) occupancy" << endl
       << "\t --subgroups      prints subgroups" << endl
       << "\t --scores         prints pwm scores" << endl
       << "\t --rate           prints rate for each gene" << endl
       << "\t --R2D            prints R for each subsequence" << endl
       << "\t --N2D            prints N for each subsequence" << endl
       << "\t --data           prints rate data for each gene" << endl
       << "\t --params         prints the parameter table" << endl
       << "\t --check-scale    checks that the scale function works with the scoring function" << endl
       << "\t --invert         if result is a data table, inverts the axes" << endl
       << "\t --gene    [name] prints only for gene with name" << endl
       << "\t --tf name [name] prints only for tf with name" << endl << endl;
  exit(1);
}
       

int main(int argc, char* argv[])
{
  int opt = 0;
  int longIndex = 0;
  string section_name("Output");
  
  string gene_name("");
  string tf_name("");
  
  bool rate       = false;
  bool sites      = false;
  bool data       = false;
  bool score      = false;
  bool maxscore   = false;
  bool scores     = false;
  bool occupancy  = false;
  bool modeocc    = false;
  bool effocc     = false;
  bool subgroups  = false;
  bool params     = false;
  bool checkscale = false;
  bool R2D        = false;
  bool N2D        = false;
  bool invert     = false;
  
  
  string infile_name;
  
  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while(opt != -1)
  {
    switch (opt)
    {
      case 'h':
        display_usage();
        break;
      case 'i':
        infile_name = optarg;
        break;
      case 's':
        section_name = optarg;
        break;
      case 0: 
        if (longOpts[longIndex].name == string("sites"))
          sites = true;
        else if (longOpts[longIndex].name == string("score"))
          score = true;
	else if (longOpts[longIndex].name == string("maxscore"))
	  maxscore = true;
        else if (longOpts[longIndex].name == string("scores"))
          scores = true;
        else if (longOpts[longIndex].name == string("rate"))
          rate = true;
        else if (longOpts[longIndex].name == string("data"))
          data = true;
        else if (longOpts[longIndex].name == string("gene"))
          gene_name = optarg;
        else if (longOpts[longIndex].name == string("tf"))
          tf_name = optarg;
        else if (longOpts[longIndex].name == string("occupancy"))
          occupancy = true;
        else if (longOpts[longIndex].name == string("modeocc"))
          modeocc = true;
        else if (longOpts[longIndex].name == string("effocc"))
          effocc = true;
        else if (longOpts[longIndex].name == string("subgroups"))
          subgroups = true;
        else if (longOpts[longIndex].name == string("section"))
          section_name = optarg;
        else if (longOpts[longIndex].name == string("params"))
          params = true;
        else if (longOpts[longIndex].name == string("check-scale"))
          checkscale = true;
        else if (longOpts[longIndex].name == string("R2D"))
          R2D = true;
        else if (longOpts[longIndex].name == string("N2D"))
          N2D = true;
        else if (longOpts[longIndex].name == string("invert"))
          invert = true;
        else
          display_usage();
        break;
      case '?':
        display_usage();
        break;
    }
    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  }
  
  if (infile_name == "")
    display_usage();
  
  ifstream infile(infile_name.c_str());

  ptree pt;
  read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);
  
  ptree& root_node   = pt.get_child("Root");
  ptree& mode_node   = root_node.get_child("Mode");
  ptree& section_node = root_node.get_child(section_name);
  
  mode_ptr mode(new Mode(infile_name,mode_node));

  mode->setVerbose(0);
  Organism embryo(section_node, mode);
  
  tfs_ptr   tfs   = embryo.getTFs();
  genes_ptr genes = embryo.getGenes();
  
  if (sites)
  {
    if      (gene_name != "" && tf_name != "")
    {
      Gene& gene = genes->getGene(gene_name);
      TF&   tf   = tfs->getTF(tf_name);
      embryo.printSites(gene, tf, cout);
    } 
    else if (gene_name != "" && tf_name == "")
    {
      Gene& gene = genes->getGene(gene_name);
      embryo.printSites(gene, cout);
    }
    else if (gene_name == "" && tf_name != "")
    {
      TF&   tf   = tfs->getTF(tf_name);
      embryo.printSites(tf, cout);
    } else
      embryo.printSites(cout);
  }
  
  if (rate)
    embryo.printRate(cout, invert);

  if (R2D)
    embryo.printR2D(cout);
  
  if (N2D)
    embryo.printN2D(cout);
  
  if (checkscale)
    embryo.checkScale(cout);

  if (data)
    embryo.printRateData(cout, invert);
  
  if (score)
    embryo.printScore(cout);

  if (maxscore)
    embryo.printMaxScore(cout);
  
  if (params)
    embryo.printParameters(cout);

  if (occupancy)
  {
    if (gene_name == "")
    {
      stringstream err;
      err << "ERROR: must specify a gene to print occupancy!" << endl;
      error(err.str());
    }
    Gene& gene = genes->getGene(gene_name);
    embryo.printOccupancy(gene, cout, invert);
  }
  
  if (scores)
  {
    if (gene_name == "")
    {
      stringstream err;
      err << "ERROR: must specify a gene to print scores!" << endl;
      error(err.str());
    }
    Gene& gene = genes->getGene(gene_name);
    embryo.printScores(gene, cout);
  }
  
  if (modeocc)
  {
    if (gene_name == "")
    {
      stringstream err;
      err << "ERROR: must specify a gene to print occupancy!" << endl;
      error(err.str());
    }
    Gene& gene = genes->getGene(gene_name);
    embryo.printModeOccupancy(gene, cout, invert);
  }
  if (effocc)
  {
    if (gene_name == "")
    {
      stringstream err;
      err << "ERROR: must specify a gene to print occupancy!" << endl;
      error(err.str());
    }
    Gene& gene = genes->getGene(gene_name);
    embryo.printEffectiveOccupancy(gene, cout, invert);
  }
  
  
  if (subgroups)
  {
    if (gene_name == "")
    {
      stringstream err;
      cerr << "ERROR: must specify a gene to print subgroups!" << endl;
      error(err.str());
    }
    Gene& gene = genes->getGene(gene_name);
    embryo.printSubgroups(gene, cout);
  }
  
  return 0;
}
