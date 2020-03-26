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

map<string,string> ligand_map()
{
  map<string,string> m;
  m["bcd"]   = "B";
  m["cad"]   = "C";
  m["dst"]   = "D";
  m["hb"]    = "H";
  m["Kr"]    = "K";
  m["kni"]   = "N";
  m["gt"]    = "G";
  m["tll"]   = "T";
  m["slp"]   = "S";
  m["prd"]   = "P";
  m["eve"]   = "E";
  m["run"]   = "R";
  m["odd"]   = "O";
  m["ftz"]   = "F";
  m["dic"]   = "I";
  m["hairy"] = "Y";
  m["moo"]   = "M";
  m["zld"]   = "Z";
  return m;
}
  
void write(ptree& out)
{
#if BOOST_VERSION / 100 % 1000 < 56
  write_xml_element(cout, 
    basic_string<ptree::key_type::value_type>(), 
    out, 
    -1, 
    boost::property_tree::xml_writer_make_settings<char>(' ', 2));
#else
  write_xml_element(cout, 
    basic_string<ptree::key_type::value_type>(), 
    out, 
    -1, 
    boost::property_tree::xml_writer_make_settings<string>(' ', 2));
#endif
}

void print_xml(Organism& organism, string& section)
{
  mode_ptr mode = organism.getMode();
  ptree out;
  ptree& system_node = out.add("System","");
  
  //ptree& description_node = system_node.add("Description", "converted from new transcpp code");
  system_node.add("Description", "converted from new transcpp code");
  
  ptree& paramset_node = system_node.add("ParamSet","");
  if (section == "Input")
    paramset_node.put("<xmlattr>.section", "input");
  else
    paramset_node.put("<xmlattr>.section", "eqparms");
  
  nuclei_ptr nuclei = organism.getNuclei();
  int nnuc = nuclei->size();
  bindings_ptr bindings = organism.getBindings();
  bindings->setMasterBindingIdx();
  
  ptree& paramlist_node = paramset_node.add("ParamList","");
  ptree& nnuc_node = paramlist_node.add("Param", "");
  nnuc_node.put("<xmlattr>.name","NoNuclei");
  nnuc_node.put("<xmlattr>.value",nnuc);
  
  ptree& af_node = paramlist_node.add("Param", "");
  af_node.put("<xmlattr>.name","BindingAffinity");
  af_node.put("<xmlattr>.value"     ,"");
  af_node.put("<xmlattr>.limit_low" ,0);
  af_node.put("<xmlattr>.limit_high",1);
  af_node.put("<xmlattr>.tweak"     ,0);
  
  ptree& fgf_node = paramlist_node.add("Param", "");
  fgf_node.put("<xmlattr>.name","Fgf");
  fgf_node.put("<xmlattr>.value"     ,1);
  fgf_node.put("<xmlattr>.limit_low" ,1);
  fgf_node.put("<xmlattr>.limit_high",1);
  fgf_node.put("<xmlattr>.tweak"     ,0);
  
  ptree& rmax_node = paramlist_node.add("Param", "");
  rmax_node.put("<xmlattr>.name","MaxRate");
  rmax_node.put("<xmlattr>.value",     255);
  rmax_node.put("<xmlattr>.limit_low" ,255);
  rmax_node.put("<xmlattr>.limit_high",255);
  rmax_node.put("<xmlattr>.tweak"     ,0);
  
  genes_ptr    genes    = organism.getGenes();
  promoter_ptr promoter = genes->getGene(0).getPromoter();
  map<string, double_param_ptr>& promoter_params = promoter->getParamMap();
  ptree& theta_node = paramlist_node.add("Param", "");
  theta_node.put("<xmlattr>.name",       "ActivationThreshold");
  theta_node.put("<xmlattr>.value",      promoter_params["Theta"]->getValue());
  theta_node.put("<xmlattr>.limit_low",  promoter_params["Theta"]->getLimLow());
  theta_node.put("<xmlattr>.limit_high", promoter_params["Theta"]->getLimHigh());
  theta_node.put("<xmlattr>.tweak",      (int) promoter_params["Theta"]->isAnnealed());
  
  ptree& q_node = paramlist_node.add("Param", "");
  q_node.put("<xmlattr>.name",       "Q");
  q_node.put("<xmlattr>.value",      promoter_params["Q"]->getValue());
  q_node.put("<xmlattr>.limit_low",  promoter_params["Q"]->getLimLow());
  q_node.put("<xmlattr>.limit_high", promoter_params["Q"]->getLimHigh());
  q_node.put("<xmlattr>.tweak",      (int) promoter_params["Q"]->isAnnealed());
  
  tfs_ptr tfs = organism.getTFs();
  int ntfs = tfs->size();
  map<string,string> m = ligand_map();
  
  ptree& ligandlist_node = paramset_node.add("LigandList","");
  
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    vector<double> coefs = tf.getCoefs();
    vector<double_param_ptr>& coef_ptrs = tf.getCoefPtrs();
    
    int ligand_type;
    if (coefs.size() == 2)
      ligand_type = 2;
    else if (coefs[0] >= 0)
      ligand_type = 0;
    else if (coefs[0] < 0)
      ligand_type = 1;
    
    ptree& ligand_node = ligandlist_node.add("Ligand","");
  
    ligand_node.put("<xmlattr>.name", tf.getName());
    ligand_node.put("<xmlattr>.ligand_id", m[tf.getName()]);
    if (coefs.size() == 2)
      ligand_node.put("<xmlattr>.ligand_type", 2);
    else if (coefs[0] >= 0)
      ligand_node.put("<xmlattr>.ligand_type", 0);
    else if (coefs[0] < 0)
      ligand_node.put("<xmlattr>.ligand_type", 1);
    ligand_node.put("<xmlattr>.include", 1);
    
    ptree& lparamlist_node = ligand_node.add("ParamList","");
    
    if (ligand_type == 0)
    {
      ptree& act_node = lparamlist_node.add("Param","");
      act_node.put("<xmlattr>.name",       "ActivationCoef");
      act_node.put("<xmlattr>.value",      coefs[0]);
      act_node.put("<xmlattr>.limit_low",  coef_ptrs[0]->getLimLow());
      act_node.put("<xmlattr>.limit_high", coef_ptrs[0]->getLimHigh());
      act_node.put("<xmlattr>.tweak",      (int) coef_ptrs[0]->isAnnealed());

    } 
    else if (ligand_type == 1)
    {
      ptree& quench_node = lparamlist_node.add("Param","");
      quench_node.put("<xmlattr>.name", "QuenchingCoef");
      quench_node.put("<xmlattr>.value",      -coefs[0]);
      quench_node.put("<xmlattr>.limit_low",  -coef_ptrs[0]->getLimLow());
      quench_node.put("<xmlattr>.limit_high", -coef_ptrs[0]->getLimHigh());
      quench_node.put("<xmlattr>.tweak",      (int) coef_ptrs[0]->isAnnealed());
      
      ptree& qd_node = lparamlist_node.add("Param","");
      qd_node.put("<xmlattr>.name",       "QuenchingDistance");
      qd_node.put("<xmlattr>.value",      100);
      qd_node.put("<xmlattr>.limit_low",  100);
      qd_node.put("<xmlattr>.limit_high", 100);
      qd_node.put("<xmlattr>.tweak",      0);
      
      ptree& qdd_node = lparamlist_node.add("Param","");
      qdd_node.put("<xmlattr>.name",       "QuenchingDistanceDelta");
      qdd_node.put("<xmlattr>.value",      50);
      qdd_node.put("<xmlattr>.limit_low",  50);
      qdd_node.put("<xmlattr>.limit_high", 50);
      qdd_node.put("<xmlattr>.tweak",      0);
      
      ptree& drep_node = lparamlist_node.add("Param","");
      drep_node.put("<xmlattr>.name",      "DirectRepressionCoef");
      drep_node.put("<xmlattr>.value",      0);
      drep_node.put("<xmlattr>.limit_low",  0);
      drep_node.put("<xmlattr>.limit_high", 0);
      drep_node.put("<xmlattr>.tweak",      0);
      
      ptree& dd_node = lparamlist_node.add("Param","");
      dd_node.put("<xmlattr>.name",       "DirectRepDistance");
      dd_node.put("<xmlattr>.value",      100);
      dd_node.put("<xmlattr>.limit_low",  100);
      dd_node.put("<xmlattr>.limit_high", 500);
      dd_node.put("<xmlattr>.tweak",      0);
      
      ptree& ddd_node = lparamlist_node.add("Param","");
      ddd_node.put("<xmlattr>.name",       "DirectRepDistanceDelta");
      ddd_node.put("<xmlattr>.value",      50);
      ddd_node.put("<xmlattr>.limit_low",  50);
      ddd_node.put("<xmlattr>.limit_high", 50);
      ddd_node.put("<xmlattr>.tweak",      0);
    }
    else if (ligand_type == 2)
    {
      ptree& quench_node = lparamlist_node.add("Param","");
      quench_node.put("<xmlattr>.name",       "QuenchingCoef");
      quench_node.put("<xmlattr>.value",      -coefs[0]);
      quench_node.put("<xmlattr>.limit_low",  -coef_ptrs[0]->getLimLow());
      quench_node.put("<xmlattr>.limit_high", -coef_ptrs[0]->getLimHigh());
      quench_node.put("<xmlattr>.tweak",      (int) coef_ptrs[0]->isAnnealed());
      
      ptree& qd_node = lparamlist_node.add("Param","");
      qd_node.put("<xmlattr>.name",       "QuenchingDistance");
      qd_node.put("<xmlattr>.value",      100);
      qd_node.put("<xmlattr>.limit_low",  100);
      qd_node.put("<xmlattr>.limit_high", 100);
      qd_node.put("<xmlattr>.tweak",      0);
      
      ptree& qdd_node = lparamlist_node.add("Param","");
      qdd_node.put("<xmlattr>.name",       "QuenchingDistanceDelta");
      qdd_node.put("<xmlattr>.value",      50);
      qdd_node.put("<xmlattr>.limit_low",  50);
      qdd_node.put("<xmlattr>.limit_high", 50);
      qdd_node.put("<xmlattr>.tweak",      0);
      
      ptree& drep_node = lparamlist_node.add("Param","");
      drep_node.put("<xmlattr>.name",      "DirectRepressionCoef");
      drep_node.put("<xmlattr>.value",     0);
      drep_node.put("<xmlattr>.limit_low",  0);
      drep_node.put("<xmlattr>.limit_high", 0);
      drep_node.put("<xmlattr>.tweak",      0);
      
      ptree& dd_node = lparamlist_node.add("Param","");
      dd_node.put("<xmlattr>.name", "DirectRepDistance");
      dd_node.put("<xmlattr>.value", 100);
      dd_node.put("<xmlattr>.limit_low",  50);
      dd_node.put("<xmlattr>.limit_high", 500);
      dd_node.put("<xmlattr>.tweak",      0);
      
      ptree& ddd_node = lparamlist_node.add("Param","");
      ddd_node.put("<xmlattr>.name",       "DirectRepDistanceDelta");
      ddd_node.put("<xmlattr>.value",      50);
      ddd_node.put("<xmlattr>.limit_low",  50);
      ddd_node.put("<xmlattr>.limit_high", 50);
      ddd_node.put("<xmlattr>.tweak",      0);
      
      ptree& act_node = lparamlist_node.add("Param","");
      act_node.put("<xmlattr>.name",       "ActivationCoef");
      act_node.put("<xmlattr>.value",      coefs[1]);
      act_node.put("<xmlattr>.limit_low",  coef_ptrs[1]->getLimLow());
      act_node.put("<xmlattr>.limit_high", coef_ptrs[1]->getLimHigh());
      act_node.put("<xmlattr>.tweak",      (int) coef_ptrs[1]->isAnnealed());
    
    }
    
    ptree& kmax_node = lparamlist_node.add("Param","");
    double_param_ptr kmax_param = tf.getKmaxParam();
    kmax_node.put("<xmlattr>.name",       "k_max");
    kmax_node.put("<xmlattr>.value",      kmax_param->getValue());
    kmax_node.put("<xmlattr>.limit_low",  kmax_param->getLimLow());
    kmax_node.put("<xmlattr>.limit_high", kmax_param->getLimHigh());
    kmax_node.put("<xmlattr>.tweak",      (int) kmax_param->isAnnealed());
    
    ptree& lambda_node = lparamlist_node.add("Param","");
    double_param_ptr lambda_param = tf.getLambdaParam();
    lambda_node.put("<xmlattr>.name",       "lambda");
    lambda_node.put("<xmlattr>.value",      lambda_param->getValue());
    lambda_node.put("<xmlattr>.limit_low",  lambda_param->getLimLow());
    lambda_node.put("<xmlattr>.limit_high", lambda_param->getLimHigh());
    lambda_node.put("<xmlattr>.tweak",      (int) lambda_param->isAnnealed());
    
    ptree& threshold_node = lparamlist_node.add("Param","");
    double_param_ptr threshold_param = tf.getThresholdParam();
    threshold_node.put("<xmlattr>.name",       "threshold");
    threshold_node.put("<xmlattr>.value",      threshold_param->getValue());
    threshold_node.put("<xmlattr>.limit_low",  threshold_param->getLimLow());
    threshold_node.put("<xmlattr>.limit_high", threshold_param->getLimHigh());
    threshold_node.put("<xmlattr>.tweak",      (int) threshold_param->isAnnealed());
    
    ptree& weightmatrix_node = ligand_node.add("WeightMatrix","");
    weightmatrix_node.put("<xmlattr>.name", tf.getName());
    weightmatrix_node.put("<xmlattr>.bsize", tf.getBindingSize());
    weightmatrix_node.put("<xmlattr>.pseudo", 1);
    weightmatrix_node.put("<xmlattr>.type", "Align");
    
    PWM& pwm = tf.getPWM();
    weightmatrix_node.put("<xmlattr>.maxscore", pwm.getMaxScore());
    weightmatrix_node.put("<xmlattr>.minscore", 0);
    vector<vector<double> > mat = pwm.getPWM(PCM);
    int pwmlen = mat.size();
    for (int j=0; j<pwmlen; j++)
    {
      ptree& position_node = weightmatrix_node.add("Position", "");
      position_node.put("<xmlattr>.Weights", to_string_(round(mat[j][0])) + ";"
                                           + to_string_(round(mat[j][1])) + ";"
                                           + to_string_(round(mat[j][2])) + ";"
                                           + to_string_(round(mat[j][3])));
    }
  }
  ptree& gintlist_node = paramset_node.add("GlobalInteractionList","");
  vector<coop_ptr>& coops = organism.getCoops()->getAllCoops();
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
  {
    ptree& coop_node = gintlist_node.add("Coop","");
    coop_node.put("<xmlattr>.actor", m[coops[i]->getTFs().first]);
    coop_node.put("<xmlattr>.target", m[coops[i]->getTFs().second]);
    coop_node.put("<xmlattr>.mode", "pair-wise");
    coop_node.put("<xmlattr>.action_type", 1);
    coop_node.put("<xmlattr>.include", 1);
    
    ptree& cparamlist = coop_node.add("ParamList","");
    ptree& kcoop = cparamlist.add("Param","");
    double_param_ptr kcoopparam = coops[i]->getKcoopParam();
    kcoop.put("<xmlattr>.name",       "Kcoop");
    kcoop.put("<xmlattr>.value",      kcoopparam->getValue());
    kcoop.put("<xmlattr>.tweak",      (int) kcoopparam->isAnnealed());
    kcoop.put("<xmlattr>.limit_high", kcoopparam->getLimHigh());
    kcoop.put("<xmlattr>.limit_low",  kcoopparam->getLimLow());
    
    map<string, double_param_ptr>& distparams = coops[i]->getDist()->getParams();
    double_param_ptr distA = distparams["A"];
    ptree& cd = cparamlist.add("Param","");
    cd.put("<xmlattr>.name",       "CoopDistance");
    cd.put("<xmlattr>.value",      distA->getValue());
    cd.put("<xmlattr>.tweak",      (int) distA->isAnnealed());
    cd.put("<xmlattr>.limit_high", distA->getLimHigh());
    cd.put("<xmlattr>.limit_low",  distA->getLimLow());
    
    ptree& cdd = cparamlist.add("Param","");
    cdd.put("<xmlattr>.name", "CoopDistanceDelta");
    cdd.put("<xmlattr>.value", 10);
    cdd.put("<xmlattr>.tweak", 0);
    cdd.put("<xmlattr>.limit_high", 10);
    cdd.put("<xmlattr>.limit_low", 10);
  }
  
  vector<coeffect_ptr>& coeffects = organism.getCoeffects()->getAllCoeffects();
  int ncoeffects = coeffects.size();
  for (int i=0; i<ncoeffects; i++)
  {
    ptree& coop_node = gintlist_node.add("Coact","");
    coop_node.put("<xmlattr>.actor", m[coeffects[i]->getActor()]);
    coop_node.put("<xmlattr>.target", m[coeffects[i]->getTarget()]);
    coop_node.put("<xmlattr>.action_type", 0);
    coop_node.put("<xmlattr>.include", 1);
    
    double_param_ptr eff = coeffects[i]->getEfficiencyParam();
    ptree& cparamlist = coop_node.add("ParamList","");
    ptree& kcoop = cparamlist.add("Param","");
    kcoop.put("<xmlattr>.name", "CoactCoef");
    kcoop.put("<xmlattr>.value", coeffects[i]->getEfficiency());
    kcoop.put("<xmlattr>.tweak", (int) eff->isAnnealed());
    kcoop.put("<xmlattr>.limit_high", eff->getLimHigh());
    kcoop.put("<xmlattr>.limit_low", eff->getLimLow());
    
    map<string, double_param_ptr>& distparams = coeffects[i]->getDist()->getParams();
    double_param_ptr distA = distparams["A"];
    double_param_ptr distB = distparams["B"];
    
    ptree& cd = cparamlist.add("Param","");
    cd.put("<xmlattr>.name", "CoactDistance");
    cd.put("<xmlattr>.value",      distA->getValue());
    cd.put("<xmlattr>.tweak",      (int) distA->isAnnealed());
    cd.put("<xmlattr>.limit_high", distA->getLimHigh());
    cd.put("<xmlattr>.limit_low",  distA->getLimLow());
    
    ptree& cdd = cparamlist.add("Param","");
    cdd.put("<xmlattr>.name", "CoactDistanceDelta");
    cdd.put("<xmlattr>.value",      distB->getValue());
    cdd.put("<xmlattr>.tweak",      (int) distB->isAnnealed());
    cdd.put("<xmlattr>.limit_high", distB->getLimHigh());
    cdd.put("<xmlattr>.limit_low",  distB->getLimLow());
  }
  
  ptree& mode_node = paramset_node.add("Mode","");
  mode_node.add("<xmlattr>.dynamic",1);
  mode_node.add("<xmlattr>.overlappingQuench",0);
  mode_node.add("<xmlattr>.overlappingCoact",0);
  mode_node.add("<xmlattr>.directRepression",0);
  mode_node.add("<xmlattr>.cooperativity",1);
  mode_node.add("<xmlattr>.scaledChisq",1);
  mode_node.add("<xmlattr>.showGroups",0);
  mode_node.add("<xmlattr>.appcomp",0);
  mode_node.add("<xmlattr>.SmoothTxnrate",1);
  mode_node.add("<xmlattr>.hilde",0);
  mode_node.add("<xmlattr>.random",0);
  
  table_ptr ratedata = organism.getRateData();
  ptree& constructlist_node = paramset_node.add("ConstructList","");
  int ngenes = genes->size();
  vector<string>& lins = ratedata->getRowNames();
  for (int i = 0; i<ngenes; i++)
  {
    ptree& construct_node = constructlist_node.add("Construct","");
    Gene& gene = genes->getGene(i);
    
    int left = gene.getLeftBound();
    construct_node.put("<xmlattr>.genotype", gene.getName());
    construct_node.put("<xmlattr>.include",  (int) gene.getInclude());
    
    scale_factor_ptr scale = gene.getScale();
    double_param_ptr scaleA = scale->getA();
    ptree& scale_node = construct_node.add("PositionEffect","");
    scale_node.put("<xmlattr>.scale_factor",scaleA->getValue());
    scale_node.put("<xmlattr>.limit_low",   scaleA->getLimLow());
    scale_node.put("<xmlattr>.limit_high",  scaleA->getLimHigh());
    scale_node.put("<xmlattr>.tweak",       (int) scaleA->isAnnealed());
    
    ptree& seq_node = construct_node.add("Sequence","");
    seq_node.put("<xmlattr>.left_bound", gene.getLeftBound());
    seq_node.put("<xmlattr>.right_bound", gene.getRightBound());
    seq_node.put("<xmlattr>.AT_bkgd", (1-mode->getGC())/2);
    seq_node.put("<xmlattr>.CG_bkgd", mode->getGC()/2);
    seq_node.put("<xmlattr>.tweak", (int) gene.getSequenceParam()->isAnnealed());
    seq_node.put("", gene.getSequenceString());
    
    ptree& sitelist_node = construct_node.add("BindingSiteList","");
    map<TF*, site_ptr_vector>& gsites = bindings->getSites(gene);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      vector<double> coefs = tf.getCoefs();
      int type = 0;
      if (coefs[0] < 0) type=1;
      if (coefs.size() == 2) type=2;
    
      site_ptr_vector& sites = gsites[&tf];
      int nsites = sites.size();
      for (int k=0; k<nsites; k++)
      {
        BindingSite& site = *(sites[k]);
        ptree& bindingsite_node = sitelist_node.add("BindingSite","");
        bindingsite_node.put("<xmlattr>.index",site.index_in_master_bindings);
        bindingsite_node.put("<xmlattr>.name",tf.getName());
        bindingsite_node.put("<xmlattr>.ligand",m[tf.getName()]);
        bindingsite_node.put("<xmlattr>.m",site.m + left);
        bindingsite_node.put("<xmlattr>.n",site.n + left);
        bindingsite_node.put("<xmlattr>.class",type);
        bindingsite_node.put("<xmlattr>.score",site.score);
        bindingsite_node.put("<xmlattr>.maxscore",tf.getMaxScore());
        bindingsite_node.put("<xmlattr>.minscore",0);
        bindingsite_node.put("<xmlattr>.k",site.K_exp_part_times_kmax);
        bindingsite_node.put("<xmlattr>.k_index",-1);
        bindingsite_node.put("<xmlattr>.k_coef",1);
        bindingsite_node.put("<xmlattr>.k_tweak",0);
        bindingsite_node.put("<xmlattr>.m_tweak",0);
        
      }
    }
   
    
    vector<double*>& col = ratedata->getCol(gene.getName());
    ptree& conclist_node = construct_node.add("MrnaConcList","");
    conclist_node.put("<xmlattr>.timepoint", 85.175000);
    int nrow = lins.size();
    for (int j=0; j<nrow; j++)
    {
      ptree& conc_node = conclist_node.add("MrnaConc","");
      conc_node.put("<xmlattr>.lin",  lins[j]);
      conc_node.put("<xmlattr>.conc", *(col[j]));
    }
      
    
    
  }
  ptree& data_node = system_node.add("Data","");
  ptree& ligandconclist_node = data_node.add("LigandConcList","");
  ligandconclist_node.put("<xmlattr>.timepoint", 85.175000);
  int nlins = lins.size();
  table_ptr tfdata = organism.getTFData();
  vector<string>& tfnames = tfdata->getColNames();
  int ntfdata = tfnames.size();
  for (int i=0; i<nlins; i++)
  {
    ptree& ligandconc_node = ligandconclist_node.add("LigandConc","");
    ligandconc_node.put("<xmlattr>.lin", lins[i]);
    for (int j=0; j<ntfdata; j++)
    {
      string name = tfnames[j];
      string id = m[name];
      ligandconc_node.put("<xmlattr>."+id, tfdata->getDataPoint("ID",lins[i],"TF",tfnames[j]));
      //cerr << "ID: " << lins[i] << endl;
      //cerr << "TF: " << tfnames[j] << endl;
      //cerr << "data: " << tfdata->getDataPoint("ID",lins[i],"TF",tfnames[j]) << endl;
    }
  }
  write(out);
}
    
    
void print_f_xml(Organism& organism, string& gname)
{
  map<string,string> lmap = ligand_map();
  ptree out;
  
  ptree& occ_node = out.add("Occupancy","");
  occ_node.put("<xmlattr>.genotype", gname);
  
  ptree& fractocclist_node = occ_node.add("FractOccList","");
  fractocclist_node.put("<xmlattr>.timepoint","85.175");
  fractocclist_node.put("<xmlattr>.lineage", "lin");
  
  ptree& bindingsitelist_node = fractocclist_node.add("BindingSiteList","");
 
  table_ptr ratedata  = organism.getRateData();
  vector<string>& ids = ratedata->getNames("ID");
  stringstream lins;
  lins << ids[0];
  int nnuc = ids.size();
  for (int i=1; i<nnuc; i++)
    lins << "," << ids[i];
  bindingsitelist_node.put("<xmlattr>.lin", lins.str());
  
  bindings_ptr bindings = organism.getBindings();
  bindings->setMasterBindingIdx();
  
  genes_ptr genes = organism.getGenes();
  tfs_ptr   tfs   = organism.getTFs();
  
  Gene& gene = genes->getGene(gname);
  map<TF*, site_ptr_vector>& gsites = bindings->getSites(gene);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    vector<double> coefs = tf.getCoefs();
    site_ptr_vector& tfsites = gsites[&tf];
    int ntfsites = tfsites.size();
    for (int j=0; j<ntfsites; j++)
    {
      BindingSite& site = *(tfsites[j]);
      ptree& bindingsite_node = bindingsitelist_node.add("BindingSite","");
      
      vector<double> fa;
      vector<double> fq;
      for (int k=0; k<nnuc; k++)
      {
        double tfa = 0;
        double tfq = 0;
        for (int l=0; l<(int) coefs.size(); l++)
        {
          if (coefs[l] > 0) tfa += site.mode_occupancy[l][k];
          if (coefs[l] < 0) tfq += site.mode_occupancy[l][k];
        }
        fa.push_back(tfa);
        fq.push_back(tfq);
      }
      stringstream ssfa;
      stringstream ssfq;
      ssfa << fa[0];
      ssfq << fq[0];
      for (int k=1; k<nnuc; k++)
      {
        ssfa << "," << fa[k];
        ssfq << "," << fq[k];
      }
      
      int type = 0;
      if (coefs[0] < 0) type=1;
      if (coefs.size() == 2) type=2;
      
      bindingsite_node.put("<xmlattr>.index", site.index_in_master_bindings);
      bindingsite_node.put("<xmlattr>.name", tf.getName());
      bindingsite_node.put("<xmlattr>.ligand", lmap[tf.getName()]);
      bindingsite_node.put("<xmlattr>.class", type);
      bindingsite_node.put("<xmlattr>.fa", ssfa.str());
      bindingsite_node.put("<xmlattr>.fq", ssfq.str());
    }
  }
  write(out);
}

void print_guts_xml(Organism& organism, string& gname)
{
  ptree out;
  table_ptr ratedata  = organism.getRateData();
  vector<string>& ids = ratedata->getNames("ID");
  int nnuc = ids.size();
  
  genes_ptr genes = organism.getGenes();
  Gene& gene = genes->getGene(gname);
  
  ptree& gutsanalysis_node = out.add("GutsAnalysis","");
  gutsanalysis_node.put("<xmlattr>.genotype", gname);
  
  ptree& gutslist_node = gutsanalysis_node.add("GutsList","");
  gutslist_node.put("<xmlattr>.timepoint", "85.175");
  gutslist_node.put("<xmlattr>.lineage", "lin");
  
  nuclei_ptr nuclei = organism.getNuclei();
  
  vector<double>& rate = nuclei->getRate(gene);
  vector<double>& N    = nuclei->getN(gene);
  
  for (int i=0; i<nnuc; i++)
  {
    ptree& guts_node = gutslist_node.add("Guts", "");
    guts_node.put("<xmlattr>.lin",ids[i]);
    guts_node.put("<xmlattr>.N",N[i]);
    guts_node.put("<xmlattr>.M",N[i]);
    guts_node.put("<xmlattr>.FGF",1);
    guts_node.put("<xmlattr>.R",rate[i]);
  }
  write(out);
}
  
// direct repression printing
void print_r_xml(Organism& organism, string& gname)
{
  map<string,string> lmap = ligand_map();
  ptree out;
  ptree& occ_node = out.add("DirectRepression","");
  occ_node.put("<xmlattr>.genotype", gname);
  
  ptree& fractocclist_node = occ_node.add("DirectRepList","");
  fractocclist_node.put("<xmlattr>.timepoint","85.175");
  fractocclist_node.put("<xmlattr>.lineage", "lin");
  
  ptree& bindingsitelist_node = fractocclist_node.add("BindingSiteList","");
 
  table_ptr ratedata  = organism.getRateData();
  vector<string>& ids = ratedata->getNames("ID");
  stringstream lins;
  lins << ids[0];
  int nnuc = ids.size();
  for (int i=1; i<nnuc; i++)
    lins << "," << ids[i];
  bindingsitelist_node.put("<xmlattr>.lin", lins.str());
  
  bindings_ptr bindings = organism.getBindings();
  bindings->setMasterBindingIdx();
  
  genes_ptr genes = organism.getGenes();
  tfs_ptr   tfs   = organism.getTFs();
  
  Gene& gene = genes->getGene(gname);
  map<TF*, site_ptr_vector>& gsites = bindings->getSites(gene);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    vector<double> coefs = tf.getCoefs();
    site_ptr_vector& tfsites = gsites[&tf];
    int ntfsites = tfsites.size();
    
    int type = 0;
    if (coefs[0] < 0) type=1;
    if (coefs.size() == 2) type=2;
    
    if (type == 0) continue; 
    
    for (int j=0; j<ntfsites; j++)
    {
      BindingSite& site = *(tfsites[j]);
      ptree& bindingsite_node = bindingsitelist_node.add("BindingSite","");
      
      vector<double> drf;
      for (int k=0; k<nnuc; k++)
      {
        double tdrf = 0;
        drf.push_back(tdrf);
      }
      stringstream ssdrf;
      ssdrf << drf[0];
      for (int k=1; k<nnuc; k++)
        ssdrf << "," << 0;
      
      bindingsite_node.put("<xmlattr>.index", site.index_in_master_bindings);
      bindingsite_node.put("<xmlattr>.name", tf.getName());
      bindingsite_node.put("<xmlattr>.ligand", lmap[tf.getName()]);
      bindingsite_node.put("<xmlattr>.class", type);
      bindingsite_node.put("<xmlattr>.drf", ssdrf.str());
    }
  }
  write(out);
}

// activation printing
void print_A_xml(Organism& organism, string& gname)
{
  map<string,string> lmap = ligand_map();
  ptree out;
  ptree& occ_node = out.add("Activation_Energy","");
  occ_node.put("<xmlattr>.genotype", gname);
  
  ptree& fractocclist_node = occ_node.add("ActEList","");
  fractocclist_node.put("<xmlattr>.timepoint","85.175");
  fractocclist_node.put("<xmlattr>.lineage", "lin");
  
  ptree& bindingsitelist_node = fractocclist_node.add("BindingSiteList","");
 
  table_ptr ratedata  = organism.getRateData();
  vector<string>& ids = ratedata->getNames("ID");
  stringstream lins;
  lins << ids[0];
  int nnuc = ids.size();
  for (int i=1; i<nnuc; i++)
    lins << "," << ids[i];
  bindingsitelist_node.put("<xmlattr>.lin", lins.str());
  
  bindings_ptr bindings = organism.getBindings();
  bindings->setMasterBindingIdx();
  
  genes_ptr genes = organism.getGenes();
  tfs_ptr   tfs   = organism.getTFs();
  
  Gene& gene = genes->getGene(gname);
  map<TF*, site_ptr_vector>& gsites = bindings->getSites(gene);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    vector<double> coefs = tf.getCoefs();
    site_ptr_vector& tfsites = gsites[&tf];
    int ntfsites = tfsites.size();
    
    int type = 0;
    if (coefs[0] < 0) type=1;
    if (coefs.size() == 2) type=2;
    
    double coef = 0;
    double midx = 0;
    if (type == 1) continue; 
    else if (type == 0) coef = coefs[0];
    else if (type == 2) { coef = coefs[1]; midx = 1; }
    
    for (int j=0; j<ntfsites; j++)
    {
      BindingSite& site = *(tfsites[j]);
      ptree& bindingsite_node = bindingsitelist_node.add("BindingSite","");
      
      vector<double> Capf;
      vector<double> Caf;
      vector<double> CaF;
      for (int k=0; k<nnuc; k++)
      {
        Capf.push_back(coef*(site.kv[k] / (1+site.kv[k])));
        Caf.push_back(coef*site.mode_occupancy[midx][k]);
        CaF.push_back(coef*site.effective_occupancy[midx][k]);
      }
      stringstream ssCapf;
      stringstream ssCaf;
      stringstream ssCaF;
      ssCapf << Capf[0];
      ssCaf  << Caf[0];
      ssCaF  << CaF[0];
      for (int k=1; k<nnuc; k++)
      {
        ssCapf << "," << Capf[k];
        ssCaf  << "," << Caf[k];
        ssCaF  << "," << CaF[k];
      }
      
      bindingsite_node.put("<xmlattr>.index", site.index_in_master_bindings);
      bindingsite_node.put("<xmlattr>.name", tf.getName());
      bindingsite_node.put("<xmlattr>.ligand", lmap[tf.getName()]);
      bindingsite_node.put("<xmlattr>.class", type);
      bindingsite_node.put("<xmlattr>.Capf", ssCapf.str());
      bindingsite_node.put("<xmlattr>.Caf", ssCaf.str());
      bindingsite_node.put("<xmlattr>.CaF", ssCaF.str());
      bindingsite_node.put("<xmlattr>.CaFFGF", ssCaF.str());
    }
  }
  write(out);
}

// effective occupancy printing
void print_F_xml(Organism& organism, string& gname)
{
  map<string,string> lmap = ligand_map();
  ptree out;
  ptree& occ_node = out.add("Occupancy","");
  occ_node.put("<xmlattr>.genotype", gname);
  
  ptree& fractocclist_node = occ_node.add("FractFOccList","");
  fractocclist_node.put("<xmlattr>.timepoint","85.175");
  fractocclist_node.put("<xmlattr>.lineage", "lin");
  
  ptree& bindingsitelist_node = fractocclist_node.add("BindingSiteList","");
 
  table_ptr ratedata  = organism.getRateData();
  vector<string>& ids = ratedata->getNames("ID");
  stringstream lins;
  lins << ids[0];
  int nnuc = ids.size();
  for (int i=1; i<nnuc; i++)
    lins << "," << ids[i];
  bindingsitelist_node.put("<xmlattr>.lin", lins.str());
  
  bindings_ptr bindings = organism.getBindings();
  bindings->setMasterBindingIdx();
  
  genes_ptr genes = organism.getGenes();
  tfs_ptr   tfs   = organism.getTFs();
  
  Gene& gene = genes->getGene(gname);
  map<TF*, site_ptr_vector>& gsites = bindings->getSites(gene);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    vector<double> coefs = tf.getCoefs();
    site_ptr_vector& tfsites = gsites[&tf];
    int ntfsites = tfsites.size();
    
    int type = 0;
    if (coefs[0] < 0) type=1;
    if (coefs.size() == 2) type=2;
    
    double midx = 0;
    if (type == 1) continue; 
    else if (type == 2) { midx = 1; }
    
    for (int j=0; j<ntfsites; j++)
    {
      BindingSite& site = *(tfsites[j]);
      ptree& bindingsite_node = bindingsitelist_node.add("BindingSite","");
      
      vector<double> F;
      for (int k=0; k<nnuc; k++)
        F.push_back(site.effective_occupancy[midx][k]);
      stringstream ssF;
      ssF << F[0];
      for (int k=1; k<nnuc; k++)
        ssF  << "," << F[k];
      
      bindingsite_node.put("<xmlattr>.index", site.index_in_master_bindings);
      bindingsite_node.put("<xmlattr>.name", tf.getName());
      bindingsite_node.put("<xmlattr>.ligand", lmap[tf.getName()]);
      bindingsite_node.put("<xmlattr>.class", type);
      bindingsite_node.put("<xmlattr>.F", ssF.str());
    }
  }
  write(out);
}

// effective occupancy printing
void print_c_xml(Organism& organism, string& gname)
{
  map<string,string> lmap = ligand_map();
  ptree out;
  ptree& occ_node = out.add("Coactivation","");
  occ_node.put("<xmlattr>.genotype", gname);
  
  ptree& fractocclist_node = occ_node.add("CoactList","");
  fractocclist_node.put("<xmlattr>.timepoint","85.175");
  fractocclist_node.put("<xmlattr>.lineage", "lin");
  
  ptree& targetlist_node = fractocclist_node.add("TargetList","");
 
  table_ptr ratedata  = organism.getRateData();
  vector<string>& ids = ratedata->getNames("ID");
  stringstream lins;
  lins << ids[0];
  int nnuc = ids.size();
  for (int i=1; i<nnuc; i++)
    lins << "," << ids[i];
  targetlist_node.put("<xmlattr>.lin", lins.str());
  
  bindings_ptr bindings = organism.getBindings();
  bindings->setMasterBindingIdx();
  
  genes_ptr genes         = organism.getGenes();
  tfs_ptr   tfs           = organism.getTFs();
  nuclei_ptr nuclei       = organism.getNuclei();
  modifying_ptr coeffects = nuclei->getCoeffects();
  Gene& gene              = genes->getGene(gname);

  map<BindingSite*, map<TF*, vector<double> > > effq = coeffects->getEffq(gene);
  map<BindingSite*, map<TF*, vector<double> > >::iterator t_it;
  map<TF*, vector<double> >::iterator                     a_it;
  
  for (t_it=effq.begin(); t_it!=effq.end(); ++t_it)
  {
    map<TF*, vector<double> >& tmap = t_it->second;
    BindingSite& target_site = *(t_it->first);
    TF&          target_tf   = *target_site.tf;
    vector<double> coefs     = target_tf.getCoefs();
    
    vector<double> fa;
    vector<double> fq;
    for (int k=0; k<nnuc; k++)
    {
      for (int l=0; l<(int) coefs.size(); l++)
      {
        double tfa = 0;
        double tfq = 0;
        if (coefs[l] > 0) tfa += target_site.mode_occupancy[l][k];
        if (coefs[l] < 0) tfq += target_site.mode_occupancy[l][k];
        fa.push_back(tfa);
        fq.push_back(tfq);
      }
    }
    stringstream ssfa;
    stringstream ssfq;
    ssfa << fa[0];
    ssfq << fq[0];
    for (int k=1; k<nnuc; k++)
    {
      ssfa << "," << fa[k];
      ssfq << "," << fq[k];
    }
    
    int type = 0;
    if (coefs[0] < 0) type=1;
    if (coefs.size() == 2) type=2;
    
    ptree& target_node = targetlist_node.add("Target","");
    target_node.put("<xmlattr>.index", target_site.index_in_master_bindings);
    target_node.put("<xmlattr>.name", target_tf.getName());
    target_node.put("<xmlattr>.ligand", lmap[target_tf.getName()]);
    target_node.put("<xmlattr>.class", type);
    target_node.put("<xmlattr>.fa", ssfa.str());
    target_node.put("<xmlattr>.fq", ssfq.str());
    
    int coac_idx = 0;
    for (a_it=tmap.begin(); a_it!=tmap.end(); ++a_it)
    {
      TF& actor = *(a_it->first);
      vector<double>& tffq = a_it->second;
      
      stringstream sstffq;
      sstffq << tffq[0];
      for (int i=1; i<(int) tffq.size(); i++)
        sstffq << "," << tffq[i];
      
      ptree& actor_node = target_node.add("Actor", "");
      
      actor_node.put("<xmlattr>.CoacIndex", coac_idx);
      actor_node.put("<xmlattr>.CoacID", lmap[actor.getName()]);
      actor_node.put("<xmlattr>.effq", sstffq.str());
    }
  }

  write(out);
}

// log odds score printing
void print_b_xml(Organism& organism, string& gname)
{
  map<string,string> lmap = ligand_map();
  ptree out;
  ptree& scores_node = out.add("PWM_Scores","");
  scores_node.put("<xmlattr>.genotype", gname);

  bindings_ptr bindings = organism.getBindings();
  genes_ptr    genes    = organism.getGenes();
  tfs_ptr      tfs      = organism.getTFs();
  Gene& gene            = genes->getGene(gname);

  stringstream ligands;
  int ntfs = tfs->size();
  ligands << tfs->getTF(0).getName();
  for (int i=1; i<ntfs; i++)
    ligands << " " << tfs->getTF(i).getName();
    
  int glength = gene.length();
  for (int i=0; i<glength; i++)
  {
    ptree& score_node = scores_node.add("Score","");
    score_node.put("<xmlattr>.index", i);
    score_node.put("<xmlattr>.pos", gene.getLeftBound() + i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      score_node.put("<xmlattr>." + tf.getName(), bindings->getScores(gene, tf).mscore[i]);
    }
  }
  scores_node.put("<xmlattr>.ligands", ligands.str());
    
  write(out);
}

void print_rate_xml(Organism& organism, string& gname)
{
  ptree out;
  table_ptr ratedata  = organism.getRateData();
  vector<string>& ids = ratedata->getNames("ID");
  int nnuc = ids.size();
  
  genes_ptr genes = organism.getGenes();
  Gene& gene = genes->getGene(gname);
  
  ptree& gutsanalysis_node = out.add("mRNARate","");
  gutsanalysis_node.put("<xmlattr>.genotype", gname);
  
  ptree& gutslist_node = gutsanalysis_node.add("MrnaConcList","");
  gutslist_node.put("<xmlattr>.timepoint", "85.175");
  gutslist_node.put("<xmlattr>.lineage", "lin");
  
  nuclei_ptr nuclei = organism.getNuclei();
  
  vector<double>& rate = nuclei->getRate(gene);
  vector<double>& N    = nuclei->getN(gene);
  
  for (int i=0; i<nnuc; i++)
  {
    ptree& guts_node = gutslist_node.add("MrnaConc", "");
    guts_node.put("<xmlattr>.lin",ids[i]);
    guts_node.put("<xmlattr>.conc",rate[i]);
  }
  write(out);
}
  
  
static const char *optString = "AbcFfhMNorsXx:";

static const struct option longOpts[] = {
    { "help",         no_argument,       NULL, 'h' },
    { "section",      required_argument, NULL, 'x' },
    { "occupancy",    no_argument,       NULL, 'f' },
    { "sequence",     no_argument,       NULL, 's' },
    { "outputxml",    no_argument,       NULL, 'o' },
    { "guts",         no_argument,       NULL, 'N' },
    { "activation",   no_argument,       NULL, 'A' },
    { "coactivation", no_argument,       NULL, 'c' },
    { "effectiveF",   no_argument,       NULL, 'F' },
    { "rate",         no_argument,       NULL, 'M' }
    //{ 0, 0, 0, 0}
};

void display_usage()
{
  cerr << endl << "\t Usage" << endl << endl
       << "\t unfold [options] input_file genotype" << endl << endl
       << "\t Options" << endl
       << "\t --help       [-h]  print this message" << endl
       << "\t --section    [-x]  use section of input file (default eqparms)" << endl
       << "\t --sequence   [-s]  use sequence-level code" << endl
       << "\t --occupancy  [-f]  prints binding site fractional occupancy" << endl
       << "\t --guts       [-N]  prints N,M,FGF,R" << endl
       << "\t --activation [-A]  prints Capf,Caf,CaF,CaFGF" << endl
       << "\t --effectiveF [-A]  prints effective occupancy F" << endl
       << "\t --logodds    [-b]  prints log odd scores for ligands" << endl
       << "\t --outputxml  [-o]  outputs the old xml format" << endl << endl;
  exit(1);
}

int main(int argc, char* argv[])
{
  int opt = 0;
  int longIndex = 0;
  string section_name("eqparms");
  
  bool sequence     = false;
  bool occupancy    = false;
  bool guts         = false;
  bool r            = false;
  bool activation   = false;
  bool coactivation = false;
  bool effectiveF   = false;
  bool logodds      = false;
  bool outputxml    = false;
  bool xmlflag      = false;
  bool rate         = false;

  string infile_name;
  
  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while(opt != -1)
  {
    switch (opt)
    {
      case 'f':
        occupancy = true;
        break;
      case 'N':
        guts = true;
        break;
      case 'A':
        activation = true;
        break;
      case 'F':
        effectiveF = true;
        break;
      case 'r':
        r = true;
        break;
      case 'c':
        coactivation = true;
        break;
      case 'h':
        display_usage();
        break;
      case 'o':
        outputxml = true;
        break;
      case 's':
        sequence = true;
        break;
      case 'b':
        logodds = true;
        break;
      case 'X':
        xmlflag = true;
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
  
  infile_name = argv[optind++];
  
  string gname("");
  bool has_gname = false;
  if (argc == (optind+1))
  {
    has_gname = true;
    gname = argv[optind];
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
  
  if (outputxml)
  {
    print_xml(embryo,section_name);
    return 0;
  }
  
  vector<string> gnames;
  if (gname == "") // doing all genes
  {
    genes_ptr genes = embryo.getGenes();
    int ngenes = genes->size();
    for (int i=0; i<ngenes; i++)
      gnames.push_back(genes->getGene(i).getName());
  }
  else
    gnames.push_back(gname);
  
  int n = gnames.size();
  for (int i=0; i<n; i++)
  {
    gname = gnames[i];
    if (occupancy)
      print_f_xml(embryo,gname);
    else if (guts)
      print_guts_xml(embryo, gname);
    else if (r)
      print_r_xml(embryo, gname);
    else if (activation)
      print_A_xml(embryo, gname);
    else if (effectiveF)
      print_F_xml(embryo, gname);
    else if (coactivation)
      print_c_xml(embryo, gname);
    else if (logodds)
      print_b_xml(embryo, gname);
    else if (rate)
      print_rate_xml(embryo, gname);
    else
      print_rate_xml(embryo, gname); // default is rate
  }
    
  return 0;
}


