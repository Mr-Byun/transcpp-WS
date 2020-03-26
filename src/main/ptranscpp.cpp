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
#include <getopt.h>
#include <libgen.h>
#include <libxml/parser.h>
#include <mpi.h>
#include "pannealer.h"
#include "move/parallelFBMove.h"
#include "intervalMix.h"
#include "unirandom.h"
#include "expHoldP.h"
#include "plsa.h"
#include "criCountP.h"
#include "dynDebug.h"

using namespace std;

using boost::property_tree::ptree;

#define to_        boost::lexical_cast
#define to_string_ boost::lexical_cast<string>

int mode_verbose;

int main(int argc, char** argv)
{
  MPIState mpiState;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiState.nnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiState.rank);
  MPI_Comm_group(MPI_COMM_WORLD, &mpiState.group);
  mpiState.comm = MPI_COMM_WORLD;

  bool isprolix = false;
  bool issteplog = false;
  bool isverbose = false;
  int iscoollog = 0;
  int readInitStates = 0;
  int optIndex;
  char *stateListFile = NULL;
  const char *readStatePrefix = NULL;
  struct option long_options[] =
    {
      {"read-state", 1, &readInitStates, 1
      },
      {"cool-log", 0, &iscoollog, 1},
      {0, 0, 0, 0}
    };
  std::string binname(basename(argv[0]));
  try
  {
    char c;
    while ((c = getopt_long(argc, argv, "plv", long_options, &optIndex)) != -1)
    {
      switch (c)
      {
      case 0:
        switch (optIndex)
        {
        case 0:
          stateListFile = optarg;
          break;
        case 1:
          break;
        default:
          throw std::runtime_error("Unrecognized option");
        }
        break;
      case 'l':
        issteplog = true;
        break;
      case 'p':
        isprolix = true;
        break;
      case 'v':
        isverbose = true;
        break;
      default:
        throw std::runtime_error("Unrecognized option");
      }
    }
    if (argc <= optind)
    {
      throw std::runtime_error("Missing input file");
    }
  }
  catch (std::exception &ex)
  {
    std::cerr << ex.what() << std::endl;
    std::cerr << "Usage: " << binname << " [ -x section_name ] input_file"
    << std::endl;
    return -1;
  }

  string xmlname(argv[optind]);
  fstream infile(xmlname.c_str());

  ptree pt;
  read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);

  ptree& root_node  = pt.get_child("Root");
  ptree& mode_node  = root_node.get_child("Mode");
  ptree& input_node = root_node.get_child("Input");
  mode_ptr mode(new Mode(xmlname, mode_node));

  mode->setVerbose(0);
  
  Organism embryo(input_node, mode);
  //embryo.printParameters(cerr);
  unsigned int seed = mode->getSeed();
  if (mode->getVerbose() >= 1)
    cerr << "Beginning annealing with seed " << seed << endl;
  unirand48 rnd;
  rnd.setSeed(seed+mpiState.rank);
  //rnd.setSeed(getpid());

  xmlDoc *doc = xmlParseFile(xmlname.c_str());
  xmlNode *docroot = xmlDocGetRootElement(doc);

  criCountP::Param frozenParam(docroot);
  intervalMix<Organism>::Param mixParam(docroot);

  pannealer<Organism, expHoldP, criCountP, parallelFBMove, intervalMix>* annealer_expHoldP;
  pannealer<Organism, plsa,     criCountP, parallelFBMove, intervalMix>* annealer_plsa;

  string bname(xmlname);
  size_t sz=bname.size();
  bname.resize(sz-4);
  string outprefix = bname + "_" +
                     ((ostringstream*)&(ostringstream() << mpiState.rank))->str();

  if (mode->getSchedule() == LAM)
  {
    //cerr << "Annealing with schedule lam" << endl;
    plsa::Param scheParam(docroot);
    annealer_plsa = new pannealer<Organism, plsa, criCountP,
                    parallelFBMove, intervalMix>(embryo, rnd, scheParam, frozenParam,
                                                 mixParam, docroot, mpiState);
    if (0 == mpiState.rank)
    {
      if (iscoollog)
        annealer_plsa->setCoolLog(file, (bname + ".log").c_str());
      if (isprolix)
        annealer_plsa->setProlix(file, (bname + ".prolix").c_str());
      if (isverbose)
      {
        annealer_plsa->setMixLog(file, (bname + ".mixlog").c_str());
      }
    }

    if (issteplog)
      annealer_plsa->setStepLog(file, (outprefix + ".steplog").c_str());
    if (readInitStates)
    {
      std::string line;
      std::ifstream is(stateListFile);
      int i = 0;
      while (!(std::getline(is,line)).eof())
      {
        if (mpiState.rank == i)
        {
          readStatePrefix = line.c_str();
          break;
        }
        ++i;
      }
      if (readStatePrefix)
        annealer_plsa->readUnifiedInitState(readStatePrefix);
      else
        throw std::runtime_error("unable to find state");
      is.close();
    }
    if (0 == mpiState.rank)
      cerr << "The energy is " << embryo.get_score() << endl;
    annealer_plsa->loop();
    xmlFreeDoc(doc);
    
    //embryo.printParameters(cerr);
    
    if (annealer_plsa->getWinner() == mpiState.rank)
    {
      cerr << "The energy is " << embryo.get_score() << " after loop" << endl;
      embryo.write("Output", root_node);
      annealer_plsa->ptreeGetResult(root_node);
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
    }
    delete annealer_plsa;
  }
  else if (mode->getSchedule() == EXP)
  {
    //cerr << "Annealing with schedule exp" << endl;
    expHoldP::Param scheParam(docroot);
    annealer_expHoldP = new pannealer<Organism, expHoldP, criCountP,
                        parallelFBMove, intervalMix>(embryo, rnd, scheParam, frozenParam,
                                                     mixParam, docroot, mpiState);
    if (0 == mpiState.rank)
    {
      if (iscoollog)
        annealer_expHoldP->setCoolLog(file, (bname + ".log").c_str());
      if (isprolix)
        annealer_expHoldP->setProlix(file, (bname + ".prolix").c_str());
      if (isverbose)
      {
        annealer_expHoldP->setMixLog(file, (bname + ".mixlog").c_str());
      }
    }

    if (issteplog)
      annealer_expHoldP->setStepLog(file, (outprefix + ".steplog").c_str());
    if (readInitStates)
    {
      std::string line;
      std::ifstream is(stateListFile);
      int i = 0;
      while (!(std::getline(is,line)).eof())
      {
        if (mpiState.rank == i)
        {
          readStatePrefix = line.c_str();
          break;
        }
        ++i;
      }
      if (readStatePrefix)
        annealer_expHoldP->readUnifiedInitState(readStatePrefix);
      else
        throw std::runtime_error("unable to find state");
      is.close();
    }
    if (0 == mpiState.rank)
      cerr << "The energy is " << embryo.get_score() << endl;
    annealer_expHoldP->loop();
    xmlFreeDoc(doc);
    
    //embryo.printParameters(cerr);
    
    if (annealer_expHoldP->getWinner() == mpiState.rank)
    {
      cerr << "The energy is " << embryo.get_score() << " after loop" << endl;
      embryo.write("Output", root_node);
      annealer_expHoldP->ptreeGetResult(root_node);
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
    }
    delete annealer_expHoldP;
  }
  
  xmlCleanupParser();
  MPI_Finalize();


  return 0;
}

