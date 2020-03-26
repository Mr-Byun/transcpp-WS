/*********************************************************************************
*                                                                                *
*     scramble.cpp                                                               *
*                                                                                *
*     Scrambles the parameters                                                   *
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

#include <unistd.h>
#include <getopt.h>
#include "unirandom.h"

using boost::property_tree::ptree;

static const char *optString = "hp::b:";

static const struct option longOpts[] = {
    { "help",        no_argument,       NULL, 'h' },
    { "permute",     optional_argument, NULL, 'p' },
    { "by",          required_argument, NULL, 'b' },
    { 0, 0, 0, 0}
};

void display_usage()
{
  cerr << endl << "\t Usage" << endl << endl
       << "\t scramble [options] template output" << endl << endl
       << "\t Options" << endl
       << "\t --help    [-h]   print this message" << endl
       << "\t --permute [-p]   permutes a datatable (defaults to RateData)" << endl
       << "\t --by      [-b]   permute the rows or cols (defaults to col)" << endl << endl;

  exit(1);
}

int main(int argc, char* argv[])
{
  int opt       = 0;
  int longIndex = 0;
  
  bool permute = false;     // by default, do not permute
  string table("RateData"); // if permuting, by default permute rate
  string by("col");         // if permuting, by default switch columns (Gene or TF, not nuc)
  
  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );

  while(opt != -1)
  {
    switch (opt)
    {
      case 'h':
        display_usage();
        break;
      case 'p':
        permute = true;
        if (optarg != NULL)
        {
          if (optarg[0] == '=') {memmove(optarg, optarg+1, strlen(optarg));}
          table = optarg;
        }
        break;
      case 'b':
        if (optarg != NULL)
        {
          if (optarg[0] == '=') {memmove(optarg, optarg+1, strlen(optarg));}
          by = optarg;
        }
        break;
      case '?':
        display_usage();
        break;
      default:
        //display_usage();
        break;
    }
    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  }
  
  string xmlname(argv[optind]);
  string xmlname2(argv[optind+1]);
  
  fstream file(xmlname.c_str());

  

  ptree pt;
  read_xml(file, pt, boost::property_tree::xml_parser::trim_whitespace);
  
  ptree& root_node  = pt.get_child("Root");
  ptree& mode_node  = root_node.get_child("Mode");
  ptree& input_node = root_node.get_child("Input");
  mode_ptr mode(new Mode(xmlname, mode_node));
  mode->setVerbose(0);

  Organism embryo(input_node, mode);
  embryo.scramble();
  if (permute)
    embryo.permute(table, by);

  //root_node.erase("Input");
  //embryo.write("Input",  root_node);
  //write_xml_element(cout, basic_string<ptree::key_type::value_type>(), output, -1, settings);
#if BOOST_VERSION / 100 % 1000 < 56
    write_xml(xmlname2, 
      pt, 
      std::locale(), 
      boost::property_tree::xml_writer_make_settings<char>(' ', 2));
#else
    write_xml(xmlname2, 
      pt, 
      std::locale(), 
      boost::property_tree::xml_writer_make_settings<string>(' ', 2));
#endif
  

  return 0;
}
