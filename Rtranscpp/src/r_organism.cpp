#include <Rcpp.h>
#include <boost/lexical_cast.hpp>
#include "r_organism.h"
#include "r_datatable.h"
#include "r_mode.h"
#include "r_parameter.h"

#define to_string_ boost::lexical_cast<string>

using namespace Rcpp;

RCPP_EXPOSED_CLASS(DataTablePtr)
RCPP_EXPOSED_CLASS(ModePtr)
RCPP_EXPOSED_CLASS(GenePtr)
RCPP_EXPOSED_CLASS(Defaults)
RCPP_EXPOSED_CLASS(TFPtr)
RCPP_EXPOSED_CLASS(DoubleParameterPtr)
RCPP_EXPOSED_CLASS(PWMParameterPtr)
RCPP_EXPOSED_CLASS(SeqParameterPtr)

OrganismPtr::OrganismPtr() :
  organism(boost::shared_ptr<Organism>(new Organism())) 
{
  mode_ptr mode(new Mode("fname"));
  mode->setVerbose(0);
  organism->setMode(mode);
}


OrganismPtr::OrganismPtr(string fname) :
  organism(boost::shared_ptr<Organism>(new Organism())) 
{
  ptree pt;
  read_xml(fname, pt, boost::property_tree::xml_parser::trim_whitespace);
  ptree& root_node  = pt.get_child("Root");
  ptree& mode_node  = root_node.get_child("Mode");
  ptree& input_node = root_node.get_child("Output");
  mode_ptr mode(new Mode(fname, mode_node));
  mode->setVerbose(0);

  organism->setMode(mode);
  organism->initialize(input_node);
}

OrganismPtr::OrganismPtr(string fname, string section) :
  organism(boost::shared_ptr<Organism>(new Organism())) 
{
  ptree pt;
  read_xml(fname, pt, boost::property_tree::xml_parser::trim_whitespace);
  ptree& root_node  = pt.get_child("Root");
  ptree& mode_node  = root_node.get_child("Mode");
  ptree& input_node = root_node.get_child(section);
  mode_ptr mode(new Mode(fname, mode_node));
  mode->setVerbose(0);
  
  organism->setMode(mode);
  organism->initialize(input_node);
}

DataTablePtr OrganismPtr::tf_data()
{
  table_ptr table = organism->getTFData();
  return DataTablePtr(table);
}

DataTablePtr OrganismPtr::rate_data()
{
  table_ptr table = organism->getRateData();
  return DataTablePtr(table);
}

ModePtr OrganismPtr::mode()
{
  mode_ptr mode = organism->getMode();
  return ModePtr(mode);
}

List OrganismPtr::get_genes()
{
  genes_ptr genes  = organism->getGenes();
  int ngenes = genes->size();
  vector<string>  gene_names(ngenes);
  
  List out(ngenes);
  for (int i=0; i<ngenes; i++)
  {
    gene_ptr gene = genes->getGeneptr(i);
    gene_names[i] = gene->getName();
    out[i] = GenePtr(gene);
  }
  out.attr("names") = gene_names;
  return out;
}

/*void OrganismPtr::set_genes(List input)
{
  genes_ptr genes;
  int ngenes = input.size();

  for (int i=0; i<ngenes; i++)
  {
    GenePtr gene = as<GenePtr>(input[i]);
    genes->add(gene.get_ptr());
  }
  organism->setGenes(genes);
}*/
  
List OrganismPtr::tfs()
{
  tfs_ptr tfs  = organism->getTFs();
  int ntfs = tfs->size();
  vector<string>  tf_names(ntfs);
  
  List out(ntfs);
  for (int i=0; i<ntfs; i++)
  {
    tf_ptr tf = tfs->getTFptr(i);
    tf_names[i] = tf->getName();
    out[i] = TFPtr(tf);
  }
  out.attr("names") = tf_names;
  return out;
}

List OrganismPtr::parameters()
{
  param_ptr_vector params = organism->getAllParameters();
  int nparams = params.size();
  List out(nparams);
  
  for (int i=0; i<nparams; i++)
  {
    iparam_ptr param = params[i];
    if (param->getType() == string("double"))
      out[i] = DoubleParameterPtr(organism,i);
    else if (param->getType() == string("PWM"))
      out[i] = PWMParameterPtr(organism,i);
    else if (param->getType() == string("Sequence"))
      out[i] = SeqParameterPtr(organism,i);
    else
      error(to_string_("No interface set for parameter of type " + param->getType()));
  }
  return out;
}

List OrganismPtr::parameter_table()
{
  param_ptr_vector params = organism->getAllParameters();
  int nparams = params.size();
  List out(4);
  CharacterVector names(nparams);
  CharacterVector types(nparams);
  CharacterVector tfs(nparams);
  LogicalVector   anneal(nparams);
  
  vector<string> col_names(4);
  vector<string> row_names(nparams);
  col_names[0] = "name";
  col_names[1] = "type";
  col_names[2] = "tf";
  col_names[3] = "anneal";
  
  for (int i=0; i<nparams; i++)
  {
    row_names[i] = to_string_(i+1);
    iparam_ptr param = params[i];
    if (param->getType() == string("double"))
    {
      DoubleParameterPtr p_ptr(organism,i);
      names[i]  = (p_ptr.name())[0];
      types[i]  = (p_ptr.type())[0];
      tfs[i]    = (p_ptr.TFName())[0];
      anneal[i] = (p_ptr.getAnneal());
      
    }
    else if (param->getType() == string("PWM"))
    {
      PWMParameterPtr p_ptr(organism, i);
      names[i]  = (p_ptr.name())[0];
      types[i]  = (p_ptr.type())[0];
      tfs[i]    = (p_ptr.TFName())[0];
      anneal[i] = (p_ptr.getAnneal());
    }
    else if (param->getType() == string("Sequence"))
    {
      SeqParameterPtr p_ptr(organism, i);
      names[i]  = (p_ptr.name())[0];
      types[i]  = (p_ptr.type())[0];
      tfs[i]    = (p_ptr.TFName())[0];
      anneal[i] = (p_ptr.getAnneal());
    }
    else
      error(to_string_("No interface set for parameter of type " + param->getType()));
  }
  out[0] = names;
  out[1] = types;
  out[2] = tfs;
  out[3] = anneal;
  
  out.attr("class")     = string("data.frame");
  out.attr("row.names") = row_names;
  out.attr("names")     = col_names;
  return out;
}

/* return the scaled rate data */
List OrganismPtr::scaled_rate_data()
{
  table_ptr         table  = organism->getRateData();
  scale_factors_ptr scales = organism->getScales();
  genes_ptr         genes  = organism->getGenes();
  
  
  vector<string>& id_names = table->getNames("ID");
  vector<string>  gene_names;
  
  int ngenes = genes->size();
  int nnuc   = id_names.size();
  
  vector<double> t_data;
  t_data.resize(nnuc);
  
  List out(ngenes);
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    const string& gname = gene.getName();
    gene_names.push_back(gname);
    scale_factor_ptr scale = gene.getScale();
    for (int j=0; j<nnuc; j++)
    {
      t_data[j] = table->getDataPoint("gene",gene.getName(),"ID",id_names[j]);
      t_data[j] = scale->scale(t_data[j]);
    }
    NumericVector new_col = wrap(t_data);
    out[i] = new_col;
  }
  out.attr("class")     = string("data.frame");
  out.attr("row.names") = id_names;
  out.attr("names")     = gene_names;
  return out;
}

/* return the calculated rates */
List OrganismPtr::rate()
{
  nuclei_ptr        nuclei = organism->getNuclei();
  genes_ptr         genes  = organism->getGenes();
  
  vector<string>& id_names = nuclei->getIDs();
  vector<string>  gene_names;
  
  int ngenes = genes->size();
  int nnuc   = id_names.size();
  
  vector<double> t_data;
  t_data.resize(nnuc);
  
  List out(ngenes);
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    const string& gname = gene.getName();
    gene_names.push_back(gname);
    vector<double>& rate = nuclei->getRate(gene);
    NumericVector new_col = wrap(rate);
    out[i] = new_col;
  }
  out.attr("class")     = string("data.frame");
  out.attr("row.names") = id_names;
  out.attr("names")     = gene_names;
  return out;
}
  

// get 2D rate data for locus with competition
List OrganismPtr::R2D()
{
  nuclei_ptr nuclei = organism->getNuclei();
  int window = nuclei->getWindow();
  int shift  = nuclei->getShift();
  
  vector<string> row_names = nuclei->getIDs();
  int nrow = nuclei->size();
  
  genes_ptr genes = nuclei->getGenes();
  int ngenes = genes->size();
  
  List out(ngenes);
  
  vector<string> gene_names;
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    int left   = gene.getLeftBound();
    gene_names.push_back(gene.getName());
    
    vector< vector<double> >& R_2D = nuclei->getR2D(gene);
    
    int nwindows = R_2D.size();
    
    vector<string> col_names;

    NumericMatrix m(nrow, nwindows);
    for (int j=0; j<nwindows; j++)
    {
      vector<double>& v = R_2D[j];
      for (int k=0; k<nrow; k++)
        m(k,j) = v[k];

      col_names.push_back(to_string_( (j+1)*shift - window/2 + left));
    }
    
    CharacterVector rows(nrow);
    CharacterVector cols(nwindows);
    for (int j=0; j<nrow; j++)
      rows[j] = row_names[j];
    for (int j=0; j<nwindows; j++)
      cols[j] = col_names[j];
    
    List dimnms = List::create(rows, cols);
    
    m.attr("dimnames") = dimnms;
    out[i] = m;
  }
  
  out.attr("names") = gene_names;
  return(out);
}

// get 2D N data for locus with competition
List OrganismPtr::N2D()
{
  nuclei_ptr nuclei = organism->getNuclei();
  int window = nuclei->getWindow();
  int shift  = nuclei->getShift();
  
  vector<string> row_names = nuclei->getIDs();
  int nrow = nuclei->size();
  
  genes_ptr genes = nuclei->getGenes();
  int ngenes = genes->size();
  
  List out(ngenes);
  
  vector<string> gene_names;
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    int left   = gene.getLeftBound();
    gene_names.push_back(gene.getName());
    
    vector< vector<double> >& N_2D = nuclei->getN2D(gene);
    
    int nwindows = N_2D.size();
    
    vector<string> col_names;

    NumericMatrix m(nrow, nwindows);
    for (int j=0; j<nwindows; j++)
    {
      vector<double>& v = N_2D[j];
      for (int k=0; k<nrow; k++)
        m(k,j) = v[k];

      col_names.push_back(to_string_( (j+1)*shift - window/2 + left));
    }
    
    CharacterVector rows(nrow);
    CharacterVector cols(nwindows);
    for (int j=0; j<nrow; j++)
      rows[j] = row_names[j];
    for (int j=0; j<nwindows; j++)
      cols[j] = col_names[j];
    
    List dimnms = List::create(rows, cols);
    
    m.attr("dimnames") = dimnms;
    out[i] = m;
  }
  
  out.attr("names") = gene_names;
  return(out);
}

// get 2D T data for locus with competition
List OrganismPtr::T2D()
{
  nuclei_ptr nuclei = organism->getNuclei();
  int window = nuclei->getWindow();
  int shift  = nuclei->getShift();
  
  vector<string> row_names = nuclei->getIDs();
  int nrow = nuclei->size();
  
  genes_ptr genes = nuclei->getGenes();
  int ngenes = genes->size();
  
  List out(ngenes);
  
  vector<string> gene_names;
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    int left   = gene.getLeftBound();
    gene_names.push_back(gene.getName());
    
    vector< vector<double> >& T_2D = nuclei->getT2D(gene);
    
    int nwindows = T_2D.size();
    
    vector<string> col_names;

    NumericMatrix m(nrow, nwindows);
    for (int j=0; j<nwindows; j++)
    {
      vector<double>& v = T_2D[j];
      for (int k=0; k<nrow; k++)
        m(k,j) = v[k];

      col_names.push_back(to_string_( (j+1)*shift - window/2 + left));
    }
    
    CharacterVector rows(nrow);
    CharacterVector cols(nwindows);
    for (int j=0; j<nrow; j++)
      rows[j] = row_names[j];
    for (int j=0; j<nwindows; j++)
      cols[j] = col_names[j];
    
    List dimnms = List::create(rows, cols);
    
    m.attr("dimnames") = dimnms;
    out[i] = m;
  }
  
  out.attr("names") = gene_names;
  return(out);
}
  
CharacterVector OrganismPtr::get_gene_names()
{
  genes_ptr genes = organism->getGenes();
  int ngenes = genes->size();
  
  CharacterVector gnames(ngenes);
  for (int i=0; i<ngenes; i++)
    gnames[i] = genes->getGene(i).getName();
  
  return gnames;
}
  
  
CharacterVector OrganismPtr::get_tf_names()
{
  tfs_ptr tfs = organism->getTFs();
  int ntfs = tfs->size();
  
  CharacterVector tfnames(ntfs);
  for (int i=0; i<ntfs; i++)
    tfnames[i] = tfs->getTF(i).getName();
  
  return tfnames;
}

CharacterVector OrganismPtr::get_nuc_names()
{
  vector<string> names = organism->getNuclei()->getIDs();
  CharacterVector out = wrap(names);
  return out;
}
  
List OrganismPtr::bindings()
{
  nuclei_ptr   nuclei   = organism->getNuclei();
  bindings_ptr bindings = nuclei->getBindings();
  genes_ptr    genes    = organism->getGenes();
  tfs_ptr      tfs      = organism->getTFs();
  
  int ngenes = genes->size();
  //int ntfs   = tfs->size();
  
  vector<string> gene_names;
  vector<string> names;
  vector<int>    tf_index;
  vector<int>    f_index;
  vector<int>    r_index;
  vector<int>    start;
  vector<int>    end;
  vector<double> score;
  vector<double> K;
  vector<double> KxKmax;
  vector<char>   strand;
    
  List out(ngenes);
  for (int i=0; i<ngenes; i++)
  {
    List gene_sites(10);
    vector<string> colnames(10);
    Gene& gene = genes->getGene(i);
    gene_names.push_back(gene.getName());
    int left_bound = gene.getLeftBound();
    vector<BindingSite*>& sites = bindings->getFsites(gene);
    
    int nsites = sites.size();
    vector<string> id_names(nsites);
    names.resize(nsites);
    tf_index.resize(nsites);
    f_index.resize(nsites);
    r_index.resize(nsites);
    start.resize(nsites);
    end.resize(nsites);
    score.resize(nsites);
    K.resize(nsites);
    KxKmax.resize(nsites);
    strand.resize(nsites);
    for (int j=0; j<nsites; j++)
    {
      id_names[j] = to_string_(j);
      BindingSite& b = *sites[j];
      names[j]    = b.tf->getName();
      tf_index[j] = b.index_in_site_map + 1;
      f_index[j]  = b.index_in_ordered_f + 1;
      r_index[j]  = b.index_in_ordered_r + 1;
      start[j]    = b.m + left_bound;
      end[j]      = b.n + left_bound;
      score[j]    = b.score;
      K[j]        = b.K_exp_part;
      KxKmax[j]   = b.K_exp_part_times_kmax;
      strand[j]   = b.orientation;
    }
    CharacterVector r_names    = wrap(names);
    IntegerVector   r_tf_index = wrap(tf_index);
    IntegerVector   r_f_index  = wrap(f_index);
    IntegerVector   r_r_index  = wrap(r_index);
    IntegerVector   r_start    = wrap(start);
    IntegerVector   r_end      = wrap(end);
    NumericVector   r_score    = wrap(score);
    NumericVector   r_K        = wrap(K);
    NumericVector   r_KxKmax   = wrap(KxKmax);
    CharacterVector r_strand   = wrap(strand);
      
    gene_sites[0] = r_names;    
    gene_sites[1] = r_tf_index;
    gene_sites[2] = r_f_index; 
    gene_sites[3] = r_r_index; 
    gene_sites[4] = r_start;   
    gene_sites[5] = r_end;     
    gene_sites[6] = r_score;   
    gene_sites[7] = r_K;       
    gene_sites[8] = r_KxKmax;  
    gene_sites[9] = r_strand;
    
    colnames[0] = "tf";
    colnames[1] = "index";
    colnames[2] = "f_index";
    colnames[3] = "r_index";
    colnames[4] = "start";
    colnames[5] = "end";
    colnames[6] = "score";
    colnames[7] = "K";
    colnames[8] = "KxKmax";
    colnames[9] = "strand";
    
    gene_sites.attr("class")     = string("data.frame");
    gene_sites.attr("row.names") = id_names;
    gene_sites.attr("names")     = colnames;
    
    out[i] = gene_sites;
  }
  out.attr("names") = gene_names;
  return(out);
}
    
List OrganismPtr::f()
{
  nuclei_ptr   nuclei   = organism->getNuclei();
  bindings_ptr bindings = nuclei->getBindings();
  genes_ptr    genes    = organism->getGenes();
  tfs_ptr      tfs      = organism->getTFs();
  
  vector<string> row_names = nuclei->getIDs();
  
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  //int nnuc   = row_names.size();
    
  vector<string> gene_names(ngenes);
  List out(ngenes);
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    gene_names[i] = gene.getName();
  
    vector<string> tf_names(ntfs);
    List gene_f(ntfs);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      tf_names[j] = tf.getName();
      site_ptr_vector sites = bindings->getSites(gene, tf);
      int nsites = sites.size();
      List tf_f(nsites);
      vector<string> site_names(nsites);
      for (int k=0; k<nsites; k++)
      {
        BindingSite& site = *(sites[k]);
        NumericVector new_col = wrap(site.total_occupancy);
        tf_f[k] = new_col;
        site_names[k] = site.tf->getName() + string("_") + to_string_(site.index_in_site_map + 1);
      }
      tf_f.attr("class") = string("data.frame");
      tf_f.attr("row.names") = row_names;
      tf_f.attr("names")     = site_names;
      gene_f[j] = tf_f;
    }
    gene_f.attr("names") = tf_names;
    out[i] = gene_f;
  }
  out.attr("names") = gene_names;
  return out;
}

List OrganismPtr::fa()
{
  nuclei_ptr   nuclei   = organism->getNuclei();
  bindings_ptr bindings = nuclei->getBindings();
  genes_ptr    genes    = organism->getGenes();
  tfs_ptr      tfs      = organism->getTFs();
  
  vector<string> row_names = nuclei->getIDs();
  
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  int ntfmodes = 0;
  for (int i=0; i<ntfs; i++)
    ntfmodes += tfs->getTF(i).getNumModes();
    
  vector<string> gene_names(ngenes);
  List out(ngenes);
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    gene_names[i] = gene.getName();
  
    vector<string> tf_names(ntfmodes);
    List gene_f(ntfmodes);
    
    int idx = 0;
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      int nmodes = tf.getNumModes();
      for (int l=0; l<nmodes; l++)
      {
        string tfname = tf.getName();
        if (l>0)
          tfname = tfname + to_string_(l);

        tf_names[idx] = tfname;
        site_ptr_vector sites = bindings->getSites(gene, tf);
        int nsites = sites.size();
        List tf_f(nsites);
        vector<string> site_names(nsites);
        for (int k=0; k<nsites; k++)
        {
          BindingSite& site = *(sites[k]);
          NumericVector new_col = wrap(site.mode_occupancy[l]);
          tf_f[k] = new_col;
          site_names[k] = site.tf->getName() + string("_") + to_string_(site.index_in_site_map + 1);
        }
        tf_f.attr("class") = string("data.frame");
        tf_f.attr("row.names") = row_names;
        tf_f.attr("names")     = site_names;
        gene_f[idx] = tf_f;
        idx++;
      }
    }
    gene_f.attr("names") = tf_names;
    out[i] = gene_f;
  }
  out.attr("names") = gene_names;
  return out;
}


List OrganismPtr::F()
{
  nuclei_ptr   nuclei   = organism->getNuclei();
  bindings_ptr bindings = nuclei->getBindings();
  genes_ptr    genes    = organism->getGenes();
  tfs_ptr      tfs      = organism->getTFs();
  
  vector<string> row_names = nuclei->getIDs();
  
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  int ntfmodes = 0;
  for (int i=0; i<ntfs; i++)
    ntfmodes += tfs->getTF(i).getNumModes();
    
  vector<string> gene_names(ngenes);
  List out(ngenes);
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    gene_names[i] = gene.getName();
  
    vector<string> tf_names(ntfmodes);
    List gene_f(ntfmodes);
    
    int idx = 0;
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      int nmodes = tf.getNumModes();
      for (int l=0; l<nmodes; l++)
      {
        string tfname = tf.getName();
        if (l>0)
          tfname = tfname + to_string_(l);

        tf_names[idx] = tfname;
        site_ptr_vector sites = bindings->getSites(gene, tf);
        int nsites = sites.size();
        List tf_f(nsites);
        vector<string> site_names(nsites);
        for (int k=0; k<nsites; k++)
        {
          BindingSite& site = *(sites[k]);
          NumericVector new_col = wrap(site.effective_occupancy[l]);
          tf_f[k] = new_col;
          site_names[k] = site.tf->getName() + string("_") + to_string_(site.index_in_site_map + 1);
        }
        tf_f.attr("class") = string("data.frame");
        tf_f.attr("row.names") = row_names;
        tf_f.attr("names")     = site_names;
        gene_f[idx] = tf_f;
        idx++;
      }
    }
    gene_f.attr("names") = tf_names;
    out[i] = gene_f;
  }
  out.attr("names") = gene_names;
  return out;
}



List OrganismPtr::scores()
{
  nuclei_ptr   nuclei   = organism->getNuclei();
  bindings_ptr bindings = nuclei->getBindings();
  genes_ptr    genes    = organism->getGenes();
  tfs_ptr      tfs      = organism->getTFs();
  
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  vector<string> gene_names(ngenes);
  vector<string> tf_names(ntfs);
  
  List out(ngenes);
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene  = genes->getGene(i);
    int glength = gene.length();
    int left    = gene.getLeftBound();
    gene_names[i] = gene.getName();
    
    vector<string> cols(ntfs);
    vector<string> rows(glength);
    for (int j=0; j<ntfs; j++)
      cols[j] = tfs->getTF(j).getName();
    for (int j=0; j<glength; j++)
      rows[j] = to_string_(j + left);
    
    List frm(3);
    vector<string> frm_names(3);
    frm_names[0] = string("forward");
    frm_names[1] = string("reverse");
    frm_names[2] = string("max");
    
    List forward(ntfs);
    List reverse(ntfs);
    List maxscore(ntfs);
    
    frm.attr("names") = frm_names;
    
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);

      TFscore& scores = bindings->getScores(gene, tf);
      NumericVector f_col = wrap(scores.fscore);
      NumericVector r_col = wrap(scores.rscore);
      NumericVector m_col = wrap(scores.mscore);
      
      forward[j]  = f_col;
      reverse[j]  = r_col;
      maxscore[j] = m_col;
    }
    forward.attr("class")  = string("data.frame");
    reverse.attr("class")  = string("data.frame");
    maxscore.attr("class") = string("data.frame");
    
    forward.attr("row.names")  = rows;
    reverse.attr("row.names")  = rows;
    maxscore.attr("row.names") = rows;
    
    forward.attr("names")  = cols;
    reverse.attr("names")  = cols;
    maxscore.attr("names") = cols;
    
    frm[0] = forward;
    frm[1] = reverse;
    frm[2] = maxscore;
    out[i] = frm;
  }
  out.attr("names") = gene_names;
  return out;
}
    

void blank() {}
  
  
RCPP_MODULE(mod_organism)
{
  class_<OrganismPtr>("Organism")
  
  .constructor()
  .constructor<string>()
  .constructor<string, string>()
  
  .property("tf_data",          &OrganismPtr::tf_data          )
  .property("rate_data",        &OrganismPtr::rate_data        )
  .property("genes",            &OrganismPtr::get_genes       , &OrganismPtr::listdummy)
  .property("tfs",              &OrganismPtr::tfs             , &OrganismPtr::listdummy)
  .property("scaled_rate_data", &OrganismPtr::scaled_rate_data )
  .property("bindings",         &OrganismPtr::bindings         )
  .property("f",                &OrganismPtr::f                )
  .property("fa",               &OrganismPtr::fa               )
  .property("F",                &OrganismPtr::F                )
  .property("mode",             &OrganismPtr::mode            , &OrganismPtr::modedummy)
  .property("R2D",              &OrganismPtr::R2D              )
  .property("N2D",              &OrganismPtr::N2D              )
  .property("T2D",              &OrganismPtr::T2D              )
  .property("rate",             &OrganismPtr::rate             )
  .property("scores",           &OrganismPtr::scores           )
  .property("score",            &OrganismPtr::score            )
  .property("par_defaults",     &OrganismPtr::get_par_defaults )
  .property("gene_names",       &OrganismPtr::get_gene_names   )
  .property("tf_names",         &OrganismPtr::get_tf_names     )
  .property("nuc_names",        &OrganismPtr::get_nuc_names    )
  .property("parameters",       &OrganismPtr::parameters      , &OrganismPtr::listdummy)
  .property("parameter_table",  &OrganismPtr::parameter_table  )
  
  .method("reset_all", &OrganismPtr::reset_all)
  .method("recalculate", &OrganismPtr::recalculate)
  ;
}
  
  
  

