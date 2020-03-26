/*********************************************************************************
*                                                                                *
*     mode.cpp                                                                   *
*                                                                                *
*     Information for how the program is run. If it is something that will       *                     
*     change the output of the program in any way, it should be here so that     *                       
*     a change in output is directly reflected in the input file. Anything       *                     
*     that does not change output (i.e. verbose level) should be handled by      *                      
*     flags, with the possible exception of modes that change how efficiently    *                        
*     the code runs.                                                             *
*                                                                                *
*********************************************************************************/

#include "mode.h"
#include "utils.h"
#include "time.h"

#include <boost/foreach.hpp>
#include <cfloat>
#include <fstream>

# define foreach_ BOOST_FOREACH



/*    Constructors    */

Mode::Mode() 
{
  scale_to         = 2000;              // the value to scale to
  min_data         = 1;                 // for functions that divide by data, the minimum to use
  p_thresh         = 0;                 // for using p-value thresholds for pwms
  gc               = 0.5;               // the default gc content to use if unspecified in input
  occupancy_method = string("dynamic"); // dynamic, combinations, sampling
  score_function   = string("sse");     // see, chisq, cc, etc
  scale_data_type  = string("area");    // what function to use when scaling data
  seed_string      = string("1000");    // the seed string in the mode
  scale_data       = true;              // do we want to scale the data?
  per_gene         = true;              // report the score per gene
  per_nuc          = true;              // report the score per nuc
  profiling        = false;             // if true, just do initial loop and exit
  self_competition = true;              // whether a TF can compete with itself
  non_specific_k   = 0;                 // adjust K for nonspecific binding energy
  verbose          = 0;                 // how much info to print during running
  num_threads      = 1;                 // if parallel, the number of threads to use
  schedule         = LAM;               // the annealing schhedule to use
  precision        = 16 ;               // the precision of printed numbers in output
  penalty_weight   = 0 ;                // by default dont use penalty
  chromatin        = false;             // dont read in chromatin by default
 
  seed = 1000; // the seed after processing
  bindingsite_list = 0;
  
  // promoter competition
  competition = false;
  window      = 500;
  shift       = 50;
  n           = 1;
  t           = 0;
}

Mode::Mode(string fname, ptree& pt) { filename = fname; read(pt); }

Mode::Mode(string fname) 
{
  filename = fname;
  fstream infile(fname.c_str());
  if (!infile.good())
    error("Could not find file with name " + fname);
  ptree pt;
  read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);
  ptree& root_node  = pt.get_child("Root");
  ptree& mode_node  = root_node.get_child("Mode");
  read(mode_node); 
}


/*    I/O   */

template<typename T>
void Mode::readNode(ptree& pt, string node_name, T* value, T def)
{
  if (pt.count(node_name))
  {
    ptree& node = pt.get_child(node_name);
    *value = node.get<T>("<xmlattr>.value", def);
  }
  else
    *value = def;
}

void Mode::readScaleData(ptree& pt)
{
  if (pt.count("ScaleData"))
  {
    ptree& scale_data_node = pt.get_child("ScaleData");
    scale_data      = scale_data_node.get<bool>(  "<xmlattr>.value", false);
    scale_data_type = scale_data_node.get<string>("<xmlattr>.type", string("area"));
    scale_to        = scale_data_node.get<double>("<xmlattr>.scale_to", 255.0);
  }
  else
    scale_data = false;
}

void Mode::process_seed_string()
{
  unsigned int out = 0;
  if (seed_string == string("filename"))
  {
    int len = filename.size();
    for (int i=0; i<len; i++)
    {
      out += (int) filename[i];
    }
  } 
  else if (seed_string == string("pid"))
    out = getpid();
  else if (seed_string == string("time"))
  {
    out = time(NULL);
  }
  else 
  {
    char *p;
    long converted = strtol(seed_string.c_str(), &p, 10);
    if (!*p)
      out = converted;
    else
      error("Seed was not set to \"filename\", \"pid\", \"time\", or a valid integer");
  }
  seed = out;
}
    
    
void Mode::readCompetition(ptree& pt)
{
  if (pt.count("Competition"))
  {
    ptree& competition_node = pt.get_child("Competition");
    competition = competition_node.get<bool>("<xmlattr>.value", false);
    window      = competition_node.get<int>("<xmlattr>.window", 500);
    shift       = competition_node.get<int>("<xmlattr>.shift", 50);
    n           = competition_node.get<int>("<xmlattr>.n", 1);
    t           = competition_node.get<double>("<xmlattr>.t", 0);
  }
  else
    competition = false;
}
    
    
  
void Mode::read(ptree& pt)
{
  ptree& mode_node = pt;
  
  readNode<int>(     mode_node, string("Verbose"),           &verbose,            0                 );
  readNode<int>(     mode_node, string("NumThreads"),        &num_threads,        1                 );
  readNode<int>(     mode_node, string("Schedule"),          &schedule,           LAM               );
  readNode<int>(     mode_node, string("Precision"),         &precision,          DBL_DIG           );
  readNode<string>(  mode_node, string("OccupancyMethod"),   &occupancy_method,   string("dynamic") );
  readNode<string>(  mode_node, string("ScoreFunction"),     &score_function,     string("sse")     );
  readNode<string>(  mode_node, string("Seed"),              &seed_string,        string("filename"));
  readNode<bool>(    mode_node, string("PerGene"),           &per_gene,           false             );
  readNode<bool>(    mode_node, string("PerNuc"),            &per_nuc,            false             );
  readNode<bool>(    mode_node, string("Profiling"),         &profiling,          false             );
  readNode<bool>(    mode_node, string("SelfCompetition"),   &self_competition,   true              );
  readNode<bool>(    mode_node, string("Chromatin"),         &chromatin,          false             );
  readNode<double>(  mode_node, string("MinData"),           &min_data,           0.0               );
  readNode<double>(  mode_node, string("PThresh"),           &p_thresh,           0.0               );
  readNode<double>(  mode_node, string("GCcontent"),         &gc,                 0.5               );
  readNode<double>(  mode_node, string("PenaltyWeight"),     &penalty_weight,     0.0               );
  readNode<double>(  mode_node, string("NonSpecificK"),      &non_specific_k,     0.0               );
  readNode<int>(     mode_node, string("BindingSiteList"),   &bindingsite_list,   0                 );
  
  readCompetition(mode_node);
  readScaleData(mode_node);
  process_seed_string();
}

void Mode::writeScaleData(ptree& pt)
{
  ptree& scale_data_node = pt.add("ScaleData","");
  scale_data_node.put("<xmlattr>.value",scale_data);
  if (scale_data)
  {
    scale_data_node.put("<xmlattr>.type",scale_data_type);
    scale_data_node.put("<xmlattr>.scale_to",scale_to);
  }
}

void Mode::write(ptree& pt)
{
  ptree& mode_node = pt.add("Mode", "");
  
  ptree& occupancy_method_node   = mode_node.add("OccupancyMethod  ", "");
  ptree& score_function_node     = mode_node.add("ScoreFunction    ", "");
  ptree& per_gene_node           = mode_node.add("PerGene          ", "");
  ptree& per_nuc_node            = mode_node.add("PerNuc           ", "");
  ptree& min_data_node           = mode_node.add("MinData          ", "");
  ptree& p_thresh_node           = mode_node.add("PThresh          ", "");
  ptree& profiling_node          = mode_node.add("Profiling        ", "");
  ptree& num_threads_node        = mode_node.add("NumThreads       ", "");
  ptree& schedule_node           = mode_node.add("Schedule         ", "");
  ptree& self_competition_node   = mode_node.add("SelfCompetition  ", "");
  ptree& precision_node          = mode_node.add("Precision        ", "");
  ptree& seed_node               = mode_node.add("Seed             ", "");
  ptree& penalty_weight_node     = mode_node.add("PenaltyWeight    ", "");
  ptree& non_specific_k_node     = mode_node.add("NonSpecificK     ", "");
  ptree& chromatin_node          = mode_node.add("Chromatin        ", "");
  ptree& bindingsite_list_node   = mode_node.add("BindingSiteList  ", "");

  occupancy_method_node.put("<xmlattr>.value", occupancy_method);
  score_function_node.put("<xmlattr>.value", score_function);
  per_gene_node.put("<xmlattr>.value", per_gene);
  per_nuc_node.put("<xmlattr>.value", per_nuc);
  min_data_node.put("<xmlattr>.value", min_data);
  p_thresh_node.put("<xmlattr>.value", p_thresh);
  profiling_node.put("<xmlattr>.value", profiling);
  num_threads_node.put("<xmlattr>.value", num_threads);
  schedule_node.put("<xmlattr>.value", schedule);
  self_competition_node.put("<xmlattr>.value", self_competition);
  precision_node.put("<xmlattr>.value", precision);
  seed_node.put("<xmlattr>.value", seed_string);
  penalty_weight_node.put("<xmlattr>.value", penalty_weight);
  non_specific_k_node.put("<xmlattr>.value", non_specific_k);
  chromatin_node.put("<xmlattr>.value", chromatin);
  bindingsite_list_node.put("<xmlattr>.value", bindingsite_list);

  writeScaleData(mode_node);
}
