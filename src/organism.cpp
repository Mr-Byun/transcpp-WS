/*********************************************************************************
*                                                                                *
*     organism.cpp                                                               *
*                                                                                *
*     An organism is defined here as a collection of nuclei. This is the master  *
*     class which holds most of the data. It creates and array of nuclei         *
*     dynamicaly based on the data present to be scored against.                 *
*                                                                                *
*********************************************************************************/


#include "organism.h"
#include <boost/foreach.hpp>
#include <limits>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/pointer_cast.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>

#define foreach_ BOOST_FOREACH
#define to_ boost::lexical_cast
#define to_string_ boost::lexical_cast<string>


/*    Constructors   */

Organism::Organism() :
    //mode(mode_ptr(new Mode)),
    distances(distances_ptr(new DistanceContainer)),
    master_tfs(tfs_ptr(new TFContainer)),
    master_genes(genes_ptr(new GeneContainer)),
    tfdata(table_ptr(new DataTable<double>)),
    ratedata(table_ptr(new DataTable<double>)),
    promoters(promoters_ptr(new PromoterContainer)),
    scale_factors(scale_factors_ptr(new ScaleFactorContainer)),
    coeffects(coeffects_ptr(new CoeffectContainer)),
    score_class(score_ptr(new Score)),
    coops(coops_ptr(new CooperativityContainer)),
    competition(competition_ptr(new Competition)),
    chromatin(chromatin_ptr(new Chromatin)),
    nuclei(nuclei_ptr(new Nuclei))
{
  move_count = 0;
  thresh = 0.5;
  test_int = 0;
}

Organism::Organism(string fname, string section) :
    //mode(mode_ptr(new Mode)),
    distances(distances_ptr(new DistanceContainer)),
    master_tfs(tfs_ptr(new TFContainer)),
    master_genes(genes_ptr(new GeneContainer)),
    tfdata(table_ptr(new DataTable<double>)),
    ratedata(table_ptr(new DataTable<double>)),
    promoters(promoters_ptr(new PromoterContainer)),
    scale_factors(scale_factors_ptr(new ScaleFactorContainer)),
    coeffects(coeffects_ptr(new CoeffectContainer)),
    score_class(score_ptr(new Score)),
    coops(coops_ptr(new CooperativityContainer)),
    competition(competition_ptr(new Competition)),
    chromatin(chromatin_ptr(new Chromatin)),
    nuclei(nuclei_ptr(new Nuclei))
{
  move_count = 0;
  thresh = 0.5;
  fstream infile(fname.c_str());
  ptree pt;
  read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);

  ptree& root_node  = pt.get_child("Root");
  ptree& mode_node  = root_node.get_child("Mode");
  ptree& input_node = root_node.get_child(section);
  mode_ptr m(new Mode(fname, mode_node));
  mode = m;
  initialize(input_node);
}

Organism::Organism(ptree& pt, mode_ptr m) :
    //mode(mode_ptr(new Mode)),
    distances(distances_ptr(new DistanceContainer)),
    master_tfs(tfs_ptr(new TFContainer)),
    master_genes(genes_ptr(new GeneContainer)),
    tfdata(table_ptr(new DataTable<double>)),
    ratedata(table_ptr(new DataTable<double>)),
    promoters(promoters_ptr(new PromoterContainer)),
    scale_factors(scale_factors_ptr(new ScaleFactorContainer)),
    coeffects(coeffects_ptr(new CoeffectContainer)),
    score_class(score_ptr(new Score)),
    coops(coops_ptr(new CooperativityContainer)),
    competition(competition_ptr(new Competition)),
    chromatin(chromatin_ptr(new Chromatin)),
    nuclei(nuclei_ptr(new Nuclei))
{
  move_count = 0;
  mode = m;
  initialize(pt);
}

void Organism::initialize(ptree& pt)
{
  move_count = 0;
  test_int = 0;
  thresh = 0.5;

  //mode->read(pt);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized the problem..." << endl;


  if (mode->getCompetition())
    competition->set(pt, mode);

  distances->setMode(mode);
  distances->add(pt);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized distance functions" << endl;

  promoters->setMode(mode);
  promoters->read(pt);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized promoter functions" << endl;

  master_tfs->add(pt, mode);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized transcription factors" << endl;

  scale_factors->setMode(mode);
  scale_factors->read(pt);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized scale factors" << endl;

  master_genes->setPromoters(promoters);
  master_genes->setScaleFactors(scale_factors);
  master_genes->read(pt);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized genes" << endl;

  ratedata->read(pt, "RateData");
  if (mode->getVerbose() >= 2)
    cerr << "Initialized rate data" << endl;

  tfdata->read(pt, "TFData");
  if (mode->getVerbose() >= 2)
    cerr << "Initialized TF data" << endl;

  coops->setMode(mode);
  coops->read(pt, distances);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized cooperativity" << endl;

  coeffects->setMode(mode);
  coeffects->read(pt,distances);
  if (mode->getVerbose() >= 2)
    cerr << "Initialized coactivation and corepression" << endl;

  chromatin->read(pt, master_genes, mode);
  
  master_tfs->setCoops(coops);
  master_tfs->setCoeffects(coeffects);

  chromatin->getParameters(params);
  competition->getParameters(params);
  distances->getParameters(params);
  master_tfs->getParameters(params);
  promoters->getParameters(params);
  scale_factors->getParameters(params);
  coops->getParameters(params);
  coeffects->getParameters(params);
  master_genes->getParameters(params);
  

  //cerr << all_params.size() << endl;
  chromatin->getAllParameters(all_params);
  //cerr << all_params.size() << endl;
  competition->getAllParameters(all_params);
  //cerr << all_params.size() << endl;
  distances->getAllParameters(all_params);
  //cerr << all_params.size() << endl;
  master_tfs->getAllParameters(all_params);
  //cerr << all_params.size() << endl;
  promoters->getAllParameters(all_params);
  //cerr << all_params.size() << endl;
  scale_factors->getAllParameters(all_params);
  //cerr << all_params.size() << endl;
  coops->getAllParameters(all_params);
  //cerr << all_params.size() << endl;
  coeffects->getAllParameters(all_params);
  //cerr << all_params.size() << endl;
  master_genes->getAllParameters(all_params);
  //cerr << all_params.size() << endl;
  


  if (mode->getVerbose() >= 2)
  {
    cerr << "Initialized parameters" << endl;
    cerr << endl;
  }

  if (mode->getVerbose() >= 1)
    printParameters(cerr);

  populate_nuclei(pt);
  if (mode->getVerbose() >= 2)
    cerr << "Created Nuclei" << endl;

  score_class->set(this);

  score();

  setMoves();
}

/*    Setters   */

void Organism::populate_nuclei(ptree& pt)
{
  nuclei->setParent(this);

  ids        = ratedata->getNames("ID");
  int nnuc   = ids.size();

  for (int i = 0; i<nnuc; i++)
  {
    string& id = ids[i];
    nuclei->addNuc(id);
  }

  nuclei->create(pt);

  if (nnuc == 0)
  {
    stringstream err;
    err << "ERROR: populate_nuclei() did not find any nuclei!" << endl;
    error(err.str());
  }
}

void Organism::populate_nuclei()
{
  nuclei->setParent(this);

  ids        = ratedata->getNames("ID");
  int nnuc   = ids.size();

  for (int i = 0; i<nnuc; i++)
  {
    string& id = ids[i];
    nuclei->addNuc(id);
  }

  nuclei->create();

  if (nnuc == 0)
  {
    stringstream err;
    err << "ERROR: populate_nuclei() did not find any nuclei!" << endl;
    error(err.str());
  }
}


double* Organism::getPrediction(Gene& gene, string& id)
{
  vector<string>& tmp_ids = nuclei->getIDs();
  int nids = tmp_ids.size();
  for (int j=0; j<nids; j++)
  {
    if (tmp_ids[j] == id)
      return &(nuclei->getRate(gene, id));
  }

  stringstream err;
  err << "ERROR: could not get prediction for " << gene.getName() << " at " << id << endl;
  error(err.str());
  return(0);
}

double* Organism::getData(Gene& gene, string& id)
{
  const string& gname = gene.getName();

  return &(ratedata->getDataPoint("gene", gname, "ID", id));
}

double* Organism::getPenalty(Gene& gene) 
{ 
  return &(nuclei->getPenalty(gene)); 
}

bindings_ptr Organism::getBindings()
{
  return nuclei->getBindings();
}




/*    Methods   */

/* scoring is actually going to be somewhat tricky in the future. Calculations
over nuclei need to be vectorized, which means that nuclei loses information about
IDs. However, at the same time we don't necessarily want nuclei and ratedata
to hold everything in the same order. This means we will need to create some map
of rate to model in scoring */

void Organism::score()
{
  score_out = score_class->getScore();
}

void Organism::printScore(ostream& os)
{
  score_class->print(os);
}

void Organism::printMaxScore(ostream& os)
{
  score_class->printMax(os);
}

void Organism::checkScale(ostream& os)
{
  score_class->checkScale(os);
}


/*    Output    */

void Organism::write(string node, ptree& pt)
{
  ptree& output = pt.add(node,"");
  if (mode->getVerbose() >= 2)
    cerr << "writing output..." << endl;
  if (mode->getCompetition())
    competition->write(output);
  distances->write(output);
  chromatin->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote distances" << endl;
  promoters->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote promoters" << endl;
  master_tfs->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote master_tfs" << endl;
  coops->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote coops" << endl;
  coeffects->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote coeffects" << endl;
  master_genes->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote master_genes" << endl;
  scale_factors->write(output);
  if (mode->getVerbose() >= 2)
    cerr << "wrote scale_factors" << endl;
  ratedata->write(output, "RateData", mode->getPrecision());
  if (mode->getVerbose() >= 2)
    cerr << "wrote ratedata" << endl;
  tfdata->write(output, "TFData", mode->getPrecision());
  if (mode->getVerbose() >= 2)
    cerr << "wrote tfdata" << endl;
  if(mode->getBindingSiteList())
    {
      // printSites(cerr);
      nuclei->write(output);
      cerr << "wrote BindingSiteList" << endl;
    }
}

void Organism::printRate(ostream& os, bool invert)
{
  int p = mode->getPrecision();   // the precision to print with
  int w = p+7; // set the minimum width, must be 6+precision to ensure columns dont merge with scientific notation
  int namew=0;
  int ngenes = master_genes->size();
  for (int i=0; i<ngenes; i++)
    namew = max(namew,(int) master_genes->getGene(i).getName().size());

  namew++;

  if (!invert)
  {
    w = max(w,namew);
    os << setw(w) << setprecision(p) << "id";
    for (int i=0; i<ngenes; i++)
      os << setw(w) << master_genes->getGene(i).getName();
    os << endl;

    vector<string>& IDs = nuclei->getIDs();
    int nids            = IDs.size();
    for (int j=0; j<nids; j++)
    {
      os << setw(w) << IDs[j];
      for (int k=0; k<ngenes; k++)
      {
        Gene& gene = master_genes->getGene(k);
        vector<double>& rate = nuclei->getRate(gene);
        os << setw(w) << rate[j];
      }
      os << endl;
    }
  }
  else
  {
    os << setw(namew) << setprecision(p) << "id";
    vector<string>& IDs = nuclei->getIDs();
    int nids            = IDs.size();
    for (int j=0; j<nids; j++)
      os << setw(w) << IDs[j];
    os << endl;

    for (int k=0; k<ngenes; k++)
    {
      Gene& gene = master_genes->getGene(k);
      os << setw(namew) << gene.getName();
      vector<string>& IDs = nuclei->getIDs();
      int nids         = IDs.size();
      vector<double>& rate = nuclei->getRate(gene);
      for (int j=0; j<nids; j++)
        os << setw(w) << rate[j];
      os << endl;
    }
  }
}


void Organism::printR2D(ostream& os)
{
  if (mode->getCompetition())
  {
    int ngenes = master_genes->size();
    for (int j=0; j<ngenes; j++)
    {
      Gene& gene = master_genes->getGene(j);
      nuclei->printR2D(gene, os);
    }
  }
  else
    error("Cannot print 2D rate with competition mode off");

}

void Organism::printN2D(ostream& os)
{
  if (mode->getCompetition())
  {
    int ngenes = master_genes->size();
    for (int j=0; j<ngenes; j++)
    {
      Gene& gene = master_genes->getGene(j);
      nuclei->printN2D(gene, os);
    }
  }
  else
    error("Cannot print 2D rate with competition mode off");
}


void Organism::printRateData(ostream& os, bool invert)
{
  int p = mode->getPrecision();   // the precision to print with
  int w = p+7; // set the minimum width, must be 6+precision to ensure columns dont merge with scientific notation
  int namew=0;
  int ngenes = master_genes->size();
  for (int i=0; i<ngenes; i++)
    namew = max(namew,(int) master_genes->getGene(i).getName().size());

  namew++;


  if (!invert)
  {
    w = max(w,namew);
    os << setw(w) << setprecision(p) << "id";
    for (int i=0; i<ngenes; i++)
    {
      Gene& gene = master_genes->getGene(i);
      if (gene.getInclude())
        os << setw(w) << gene.getName();
    }
    os << endl;

    vector<string>& IDs = nuclei->getIDs();
    int nids         = IDs.size();
    for (int j=0; j<nids; j++)
    {
      string& id = IDs[j];
      os << setw(w) << id;
      for (int k=0; k<ngenes; k++)
      {
        Gene& gene = master_genes->getGene(k);
        if (!gene.getInclude()) continue;
        scale_factor_ptr scale = gene.getScale();
        double datapoint = ratedata->getDataPoint("gene",gene.getName(),"ID",id);
        datapoint = scale->scale(datapoint);
        os << setw(w) << datapoint;
      }
      os << endl;
    }
  }
  else
  {
    os << setw(namew) << setprecision(p) << "id";

    vector<string>& IDs = nuclei->getIDs();
    int nids         = IDs.size();
    for (int j=0; j<nids; j++)
      os << setw(w) << IDs[j];

    os << endl;

    for (int k=0; k<ngenes; k++)
    {
      Gene& gene = master_genes->getGene(k);
      if (!gene.getInclude()) continue;
      scale_factor_ptr scale = gene.getScale();
      os << setw(namew) << gene.getName();

      vector<string>& IDs = nuclei->getIDs();
      int nids         = IDs.size();
      for (int j=0; j<nids; j++)
      {
        string& id = IDs[j];
        double datapoint = ratedata->getDataPoint("gene",gene.getName(),"ID",id);
        datapoint = scale->scale(datapoint);
        os << setw(w) << datapoint;
      }
      os << endl;
    }
  }
}

void Organism::printSites(ostream& os)
{
  nuclei->printSites(os);
}

void Organism::printSites(Gene& gene, ostream& os)
{
  nuclei->printSites(gene, os);
}

void Organism::printSites(TF& tf, ostream& os)
{
  nuclei->printSites(tf, os);
}

void Organism::printSites(Gene& gene, TF& tf, ostream& os)
{
  nuclei->printSites(gene, tf, os);
}

void Organism::printScores(Gene& gene, ostream& os)
{
  nuclei->printScores(gene, os);
}

void Organism::printSubgroups(Gene& gene, ostream& os)
{
  nuclei->printSubgroups(gene, os);
}

void Organism::printOccupancy(Gene& gene, ostream& os, bool invert)
{
  nuclei->printOccupancy(gene, os, invert);
}

void Organism::printModeOccupancy(Gene& gene, ostream& os, bool invert)
{
  nuclei->printModeOccupancy(gene, os, invert);
}

void Organism::printEffectiveOccupancy(Gene& gene, ostream& os, bool invert)
{
  nuclei->printEffectiveOccupancy(gene, os, invert);
}

/*    Move    */

/* scramble now changes the input node as well as the parameter value, so that
the node can just be reprinted without reprinting anything else */

void Organism::scramble()
{
  boost::minstd_rand baseGen(getpid());
  boost::uniform_real<> uniDblUnit(0,1);
  boost::variate_generator<boost::minstd_rand&, boost::uniform_real<> > uniDblGen(baseGen, uniDblUnit);

  uniDblGen();
  int nparams = params.size();
  for (int i=0; i<nparams; i++)
    params[i]->scramble(uniDblGen());
  
  ResetAll();
  score();
}

void Organism::permute(string& table, string& by)
{
  if (table == string("RateData"))
    ratedata->permute(by, mode->getPrecision());
  else if (table == string("TFData"))
    tfdata->permute(by, mode->getPrecision());
  else
  {
    stringstream err;
    err << "ERROR: no data table with name " << table << endl;
    error(err.str());
  }
}

/***************  Move Generation  **********************************************/

void Organism::setMoves()
{
  moves.clear();
  restores.clear();
  all_moves.clear();
  all_restores.clear();

  setPVectorMoves(moves, restores, params);
  setPVectorMoves(all_moves, all_restores, all_params);
  
}
  
void Organism::setPVectorMoves(vector<boost::function<void (Organism*)> >& mvec, vector<boost::function<void (Organism*)> >& rvec, param_ptr_vector& pvec)
{
  string move;
  int npvec = pvec.size();
  for (int i=0; i<npvec; i++)
  {
    move = pvec[i]->getMove();

    if (move == string("Scores"))
    {
      TF& tf = master_tfs->getTF(pvec[i]->getTFName());
      mvec.push_back(boost::bind(&Organism::moveScores, this, boost::ref(tf)));
      rvec.push_back(boost::bind(&Organism::restoreScores, this, boost::ref(tf)));
    }
    else if (move == string("PWM"))
    {
      TF& tf = master_tfs->getTF(pvec[i]->getTFName());
      mvec.push_back(boost::bind(&Organism::movePWM, this, boost::ref(tf)));
      rvec.push_back(boost::bind(&Organism::restorePWM, this, boost::ref(tf)));
    }
    else if (move == string("Sites"))
    {
      TF& tf = master_tfs->getTF(pvec[i]->getTFName());
      mvec.push_back(boost::bind(&Organism::moveSites, this, boost::ref(tf)));
      rvec.push_back(boost::bind(&Organism::restoreSites, this, boost::ref(tf)));
    }
    else if (move == string("Lambda"))
    {
      TF& tf = master_tfs->getTF(pvec[i]->getTFName());
      mvec.push_back(boost::bind(&Organism::moveLambda, this, boost::ref(tf)));
      rvec.push_back(boost::bind(&Organism::restoreLambda, this, boost::ref(tf)));
    }
    else if (move == string("Kmax"))
    {
      TF& tf = master_tfs->getTF(pvec[i]->getTFName());
      mvec.push_back(boost::bind(&Organism::moveKmax, this, boost::ref(tf)));
      rvec.push_back(boost::bind(&Organism::restoreKmax, this, boost::ref(tf)));
    }
    else if (move == string("Coef"))
    {
      double_param_ptr p = boost::dynamic_pointer_cast<Parameter<double> >(pvec[i]);
      mvec.push_back(boost::bind(&Organism::moveCoef, this, boost::ref(p)));
      rvec.push_back(boost::bind(&Organism::restoreCoef, this, boost::ref(p)));
    }
    else if (move == string("CoopD"))
    {
      mvec.push_back(boost::bind(&Organism::moveCoopD, this));
      rvec.push_back(boost::bind(&Organism::restoreCoopD, this));
    }
    else if (move == string("Kcoop"))
    {
      mvec.push_back(boost::bind(&Organism::moveKcoop, this));
      rvec.push_back(boost::bind(&Organism::restoreKcoop, this));
    }
    else if (move == string("Quenching"))
    {
      mvec.push_back(boost::bind(&Organism::moveQuenching, this));
      rvec.push_back(boost::bind(&Organism::restoreQuenching, this));
    }
    else if (move == string("QuenchingCoef"))
    {
      mvec.push_back(boost::bind(&Organism::moveQuenchingCoef, this));
      rvec.push_back(boost::bind(&Organism::restoreQuenchingCoef, this));
    }
    else if (move == string("Coeffect"))
    {
      mvec.push_back(boost::bind(&Organism::moveCoeffect, this));
      rvec.push_back(boost::bind(&Organism::restoreCoeffect, this));
    }
    else if (move == string("CoeffectEff"))
    {
      mvec.push_back(boost::bind(&Organism::moveCoeffectEff, this));
      rvec.push_back(boost::bind(&Organism::restoreCoeffectEff, this));
    }
    else if (move == string("Window"))
    {
      mvec.push_back(boost::bind(&Organism::moveWindow, this));
      rvec.push_back(boost::bind(&Organism::moveWindow, this));
    }
    else if (move == string("Promoter"))
    {
      mvec.push_back(boost::bind(&Organism::movePromoter, this));
      rvec.push_back(boost::bind(&Organism::movePromoter, this));
    }
    else if (move == string("Kacc"))
    {
      mvec.push_back(boost::bind(&Organism::moveKacc, this));
      rvec.push_back(boost::bind(&Organism::restoreKacc, this));
    }
    else if (move == string("Null"))
    {
      mvec.push_back(boost::bind(&Organism::null_function, this));
      rvec.push_back(boost::bind(&Organism::null_function, this));
    }
    else if (move == string("ResetAll"))
    {
      if (pvec[i]->isAnnealed())
        warning("Move function for parameter " + pvec[i]->getParamName() + " is ResetAll. Annealing may be very slow");
      mvec.push_back(boost::bind(&Organism::ResetAll, this));
      rvec.push_back(boost::bind(&Organism::ResetAll, this));
    }
    else
      error("setMvec() Could not find move function with name " + move + " for parameter " + pvec[i]->getParamName());
  }
}


/*
Now we define the move and restore functions. I have tried to order these
according to how much is reset in each case. I thought this might be easier to
maintain if I set a "reset level" for each parameter, however there are no clear
rules for what needs to be reset depending on what parameter we tweak. Instead,
this requires pretty good knowledge of the transcription equations themselves, 
so be careful when adding to these. 
 
There are four types of functions:
 
save     - saves old data structures
update   - implements any preprocessing of interactions and memory managements
calc     - preprocessing is unnecessary, so simply rerun the calculation
restore  - copies over old data structures 
*/

/* The slowest way to reset everything, this function deletes all calculated objects
and repopulates them from scratch! */
void Organism::Recalculate()
{
  nuclei = nuclei_ptr(new Nuclei);
  
  master_tfs->setCoops(coops);
  master_tfs->setCoeffects(coeffects);

  params.clear();
  all_params.clear();
  
  chromatin->getParameters(params);
  competition->getParameters(params);
  distances->getParameters(params);
  master_tfs->getParameters(params);
  promoters->getParameters(params);
  scale_factors->getParameters(params);
  coops->getParameters(params);
  coeffects->getParameters(params);
  master_genes->getParameters(params);
  

  chromatin->getAllParameters(all_params);
  competition->getAllParameters(all_params);
  distances->getAllParameters(all_params);
  master_tfs->getAllParameters(all_params);
  promoters->getAllParameters(all_params);
  scale_factors->getAllParameters(all_params);
  coops->getAllParameters(all_params);
  coeffects->getAllParameters(all_params);
  master_genes->getAllParameters(all_params);
  
  populate_nuclei();
  
  score_class->set(this);

  score();

  setMoves();
}

/* clear all data and reset from scratch. If you suspect a move function
is not working, this is the easiest way to check, however it should generally
be avoided as it will be very slow */

void Organism::ResetAll()
{
  if (mode->getVerbose() >= 3)
    cerr << "Reseting everything" << endl;
  
  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->updateScores(gene);
    nuclei->updateSites(gene);
    nuclei->updateSubgroups(gene);
    nuclei->updateCoeffects(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    if (mode->getCompetition() == true)
      nuclei->resizeWindow(gene);
    nuclei->updateR(gene);

  }
}

/* this will be necessary if something changes the way we score sequence, for instance
tweaking pwms. In general we believe this is a bad idea, but it could be useful
for comparing to other results (segal) or in very careful applications. If you
decide to use this feature, be ready to defend your decision! */
void Organism::moveScores(TF& tf)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving sequence scores" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->saveSites(gene,tf);
    nuclei->saveScores(gene,tf);
    nuclei->saveSubgroups(gene);
    nuclei->saveCoeffects(gene);
    nuclei->saveQuenching(gene);

    nuclei->updateScores(gene,tf);
    nuclei->updateSites(gene,tf);
    nuclei->updateSubgroups(gene);
    nuclei->updateCoeffects(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::restoreScores(TF& tf)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring sequence scores" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreScores(gene,tf);
    nuclei->restoreSites(gene,tf);
    nuclei->restoreAllOccupancy(gene);
    nuclei->restoreSubgroups(gene);
    nuclei->restoreCoeffects(gene);
    nuclei->restoreQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

/* this will be necessary if something changes the way we score sequence, for instance
tweaking pwms. In general we believe this is a bad idea, but it could be useful
for comparing to other results (segal) or in very careful applications. If you
decide to use this feature, be ready to defend your decision! */
void Organism::movePWM(TF& tf)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving PWM" << endl;

  int ngenes  = master_genes->size();

  //params[idx]->print(cerr);
  
#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->saveSites(gene, tf);
    nuclei->saveScores(gene, tf);
    nuclei->saveSubgroups(gene);
    nuclei->saveCoeffects(gene);
    nuclei->saveQuenching(gene);

    nuclei->updateScores(gene,tf);
    nuclei->updateSites(gene,tf);
    nuclei->updateSubgroups(gene);
    nuclei->updateCoeffects(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

void Organism::restorePWM(TF& tf)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring PWM" << endl;
  
  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreScores(gene,tf);
    nuclei->restoreSites(gene,tf);
    nuclei->restoreAllOccupancy(gene);
    nuclei->restoreSubgroups(gene);
    nuclei->restoreCoeffects(gene);
    nuclei->restoreQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}


/* if thresholds are changed, the sites will need to be repopulated for the tf
changed */
void Organism::moveSites(TF& tf)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving Sites" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->saveSites(gene, tf);
    nuclei->saveSubgroups(gene);
    nuclei->saveCoeffects(gene);
    nuclei->saveQuenching(gene);

    nuclei->updateSites(gene,tf);
    nuclei->updateSubgroups(gene);
    nuclei->updateCoeffects(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

void Organism::restoreSites(TF& tf)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring Sites" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreSites(gene, tf);
    //nuclei->updateSites(gene);
    nuclei->restoreAllOccupancy(gene);
    nuclei->restoreSubgroups(gene);
    nuclei->restoreCoeffects(gene);
    nuclei->restoreQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

/* if we move cooperativity distance we need to repopulate subgroups, but
quenching and coeffect interactions are unchanged */
void Organism::moveCoopD()
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving cooperativity distance" << endl;


  int ngenes  = master_genes->size();

#ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->saveSubgroups(gene);

    nuclei->updateSubgroups(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

void Organism::restoreCoopD()
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring cooperativity distance" << endl;


  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreAllOccupancy(gene);
    nuclei->restoreSubgroups(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

/* If we change lambda we need dont need to scan for new sites, but we do
need to update K and recalc occupancy and interations */
void Organism::moveLambda(TF& tf)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving lambda" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    // dont bother saving K and Lambda since that takes nearly as long as calculating
    nuclei->updateKandLambda(gene, tf);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

void Organism::restoreLambda(TF& tf)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring lambda" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreAllOccupancy(gene);
    nuclei->updateKandLambda(gene, tf);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::moveKacc()
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving kacc" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    // dont bother saving K and Lambda since that takes nearly as long as calculating
    nuclei->updateKandLambda(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

void Organism::restoreKacc()
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring kacc" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreAllOccupancy(gene);
    nuclei->updateKandLambda(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

// Moving Kmax is about the same as moving lambda. Maybe I should just ignore..
void Organism::moveKmax(TF& tf)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving kmax for tf " << tf.getName() << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->updateK(gene, tf);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

void Organism::restoreKmax(TF& tf)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring kmax" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreAllOccupancy(gene);
    nuclei->updateK(gene, tf);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }

}

// if we move cooperativity we simply need to redo occupancy calculations
void Organism::moveKcoop()
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving kcoop" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::restoreKcoop()
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring kcoop" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreAllOccupancy(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

/* If we change coeffect distances we need to repopulate the coeffects*/
void Organism::moveCoeffect()
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving coeffects" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveModeOccupancy(gene);
    nuclei->saveQuenching(gene);
    nuclei->saveCoeffects(gene);

    nuclei->updateCoeffects(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::restoreCoeffect()
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring coeffects" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreCoeffects(gene);
    nuclei->restoreQuenching(gene);
    nuclei->restoreModeOccupancy(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

/* If we allow a coefficient to switch from activator to repressor, we use this
function. If the bounds dont allow this, you can simply point the move generator
to move quenching or move activations */
void Organism::moveCoef(double_param_ptr p)
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving TF coefficient" << endl;
  
  val  = p->getValue();
  prev = p->getPrevious();
  
  if ( val >= 0 && prev >= 0)
    movePromoter();
  else if ( val <= 0 && prev <= 0)
    moveQuenchingCoef();
  else
    moveQuenching();
}

void Organism::restoreCoef(double_param_ptr p)
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring TF coefficient" << endl;
  
  if ( val >= 0 && prev >= 0)
    movePromoter();
  else if ( val <= 0 && prev <= 0)
    restoreQuenchingCoef();
  else
    restoreQuenching();
}

/* These functions are used if a quencher has been added or removed */
void Organism::moveQuenching()
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving quenching" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveAllOccupancy(gene);
    //nuclei->saveSites(gene, tf);
    //nuclei->saveScores(gene, tf);
    //nuclei->saveSubgroups(gene);
    nuclei->saveCoeffects(gene);
    nuclei->saveQuenching(gene);

    //nuclei->updateScores(gene,tf);
    //nuclei->updateSites(gene,tf);
    //nuclei->updateSubgroups(gene);
    nuclei->updateCoeffects(gene);
    nuclei->updateQuenching(gene);
    nuclei->calcOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    nuclei->updateR(gene);
  }
}

void Organism::restoreQuenching()
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring quenching" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif
    
  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    //nuclei->restoreScores(gene,tf);
    //nuclei->restoreSites(gene,tf);
    nuclei->restoreAllOccupancy(gene);
    //nuclei->restoreSubgroups(gene);
    nuclei->restoreCoeffects(gene);
    nuclei->restoreQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}


/* These functions are used if a the number of quenchers is unchanged */
void Organism::moveQuenchingCoef()
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving quenching coefficient" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveEffectiveOccupancy(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::restoreQuenchingCoef()
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring quenching coefficient" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreEffectiveOccupancy(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

/* If we have moved coactivation or corepression efficiency */
void Organism::moveCoeffectEff()
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving coeffect coefficient" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->saveModeOccupancy(gene);
    nuclei->calcCoeffects(gene);
    nuclei->calcQuenching(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

void Organism::restoreCoeffectEff()
{
  if (mode->getVerbose() >= 3)
    cerr << "Restoring coeffect coefficient" << endl;

  int ngenes  = master_genes->size();

#ifdef PARALLEL
    #pragma omp parallel for num_threads(mode->getNumThreads())
    #endif

  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->restoreModeOccupancy(gene);
    //nuclei->updateN();
    nuclei->updateR(gene);
  }
}

/* if we have moved the promoter properties we simply call this */
void Organism::movePromoter()
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving promoter parameter" << endl;

  int ngenes  = master_genes->size();

  #ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif
  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->updateR(gene);
  }
}

/* if we have moved the promoter properties we simply call this */
void Organism::moveWindow()
{
  if (mode->getVerbose() >= 3)
    cerr << "Moving window parameter" << endl;

  int ngenes  = master_genes->size();

  #ifdef PARALLEL
  #pragma omp parallel for num_threads(mode->getNumThreads())
  #endif
  for (int j=0; j<ngenes; j++)
  {
    Gene& gene = master_genes->getGene(j);
    if (!gene.getInclude()) continue;
    nuclei->resizeWindow(gene);
    nuclei->updateR(gene);
  }
}

/* this is somewhat awkward, but we dont need to do anything but rescore if
scale factors are used, so this is a null funtions */

void Organism::null_function()
{}


/*    Annealing   */

int Organism::getStateSize()
{
  int nparams = params.size();
  int size = 0;
  for (int i=0; i<nparams; i++)
    size += params[i]->getSize();
  
  return size;
}

void Organism::serialize(void * buf) const
{
  char * cbuf = (char *) buf;
  int nparams = params.size();
  for (int i = 0; i < nparams; ++i)
  {
    params[i]->serialize(cbuf);
    cbuf += params[i]->getSize();
  }   
}

void Organism::deserialize(void const *buf)
{
  char const * cbuf = (char const *) buf;
  int nparams = params.size();
  for (int i = 0; i < nparams; ++i)
  {
    params[i]->deserialize(cbuf);
    cbuf += params[i]->getSize();
  }
  ResetAll();
  score();
}



int    Organism::getDimension() const
{
  return params.size();
}

double Organism::get_score()
{
  return score_out;
}

void   Organism::generateMove(int idx, double theta)
{
  move_count++;
  /*if (move_count == 1)
    adjustThresholds(thresh);
  if (move_count > 50000)
  {
    if (move_count % 5000 == 0)
    {
      if (thresh != 0.0)
      {
        thresh = max(0.0, thresh-0.005);
        adjustThresholds(thresh);
      }
    }
  }*/
  
  //if (move_count % 1000 == 0)
  //  cerr << "moves: " << move_count << "  thresh: " << thresh << endl;
  
  previous_score_out = score_out;

  if (mode->getVerbose() >= 3)
    cerr << "generating move..." << endl;

  params[idx]->tweak(theta);

  if (!params[idx]->isOutOfBounds())
  {
    moves[idx](this);
    score();
  }
  else
  {
    if (mode->getVerbose() >= 3)
      cerr << "move is out of bounds!" << endl;
    score_out = numeric_limits<double>::max();
    params[idx]->restore();
  }


}

void Organism::adjustThresholds(double percent)
{
  cerr << "adjusting thresholds to " << percent*100 << " percent" << endl;
  int ntfs = master_tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = master_tfs->getTF(i);
    double max_score = tf.getMaxScore();
    double new_threshold = max_score * percent;
    tf.setThreshold(new_threshold);
    //cerr << "threshold for " << tf.getName() << " is " << new_threshold << endl;
  }
  ResetAll();
}
  
void Organism::move(int idx)
{
  moves[idx](this);
  score();
}

void Organism::move_all(int idx)
{
  
  all_moves[idx](this);
  score();
}

void Organism::restore_all(int idx)
{
  all_restores[idx](this);
  score();
}

void   Organism::restoreMove(int idx)
{
  if (mode->getVerbose() >= 3)
    cerr << "restoring move" << endl;

  params[idx]->restore();

  // we only want to restore a move if it was out of bounds
  if (!params[idx]->isOutOfBounds())
  {
    restores[idx](this);
  }

  score();
  if (mode->getVerbose() >= 3)
    cerr << "the score is now: " << setprecision(16) << score_out << endl;
}

void Organism::printParameters(ostream& os)
{
  os << "Parameters" << endl;
  int nparams = params.size();

  for (int i=0; i<nparams; i++)
  {
    if (i==0)
    {
      os << setw(5) << "idx";
      params[i]->printHeader(os);
    }
    os << setw(5) << i;
    params[i]->print(os);
  }
  os << endl;
}

string Organism::getParamName(int idx)
{
  stringstream out;

  if (params[idx]->is_tf_param())
    out << params[idx]->getTFName() << "_";
  out << params[idx]->getParamName();

  return out.str();

}

vector<double>& Organism::getN(int gidx)
{
  Gene& gene = master_genes->getGene(gidx);
  vector<double>& out = nuclei->getN(gene);
  return out;
}

vector<double>& Organism::getR(int gidx)
{
  Gene& gene = master_genes->getGene(gidx);
  vector<double>& out = nuclei->getRate(gene);
  return out;
}
