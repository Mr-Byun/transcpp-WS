/*********************************************************************************
*                                                                                *
*     bindings.cpp                                                               *
*                                                                                *
*     Constains a map of all site on all genes                                   *
*                                                                                *
*********************************************************************************/

#include "bindings.h"

#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/lexical_cast.hpp>
#include <limits>

#define foreach_ BOOST_FOREACH
#define to_string_ boost::lexical_cast<string>

Bindings::Bindings() {nnuc=0;}

void Bindings::clear()
{
  scores.clear();
  saved_scores.clear();
  sites.clear();
  saved_sites.clear();
}

void Bindings::setGenes(genes_ptr g)
{
  genes = g;
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    gene_scores_map_ptr gscores(new gene_scores_map);
    gene_scores_map_ptr saved_gscores(new gene_scores_map);
    
    gene_sites_map_ptr gsites(new gene_sites_map);
    gene_sites_map_ptr saved_gsites(new gene_sites_map);
    
    scores[&gene]       = gscores;
    saved_scores[&gene] = saved_gscores;
    
    sites[&gene]       = gsites;
    saved_sites[&gene] = saved_gsites;
  }
}
    
map<TF*, site_ptr_vector>& Bindings::getSites(Gene& gene)
{
  return *(sites[&gene]);
}

site_ptr_vector& Bindings::getSites(Gene& gene, TF& tf) 
{
  gene_sites_map& gsites = *(sites[&gene]);
  return gsites[&tf];
}

void Bindings::setMasterBindingIdx()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  int idx = 0;
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    gene_sites_map& gsites = *(sites[&gene]);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      site_ptr_vector& tfsites = gsites[&tf];
      int nsites = tfsites.size();
      for (int k=0; k<nsites; k++)
        tfsites[k]->index_in_master_bindings = idx++;
    }
  }
}
        
      
  
/* these functions are broken, but is also never used, so commented out for now
bool Bindings::hasScores(Gene& gene, TF& tf)
{
  if (scores.find(&gene) == scores.end())
    return false;
  else
  {
    if (scores[&gene].find(&tf) == scores[&gene].end())
      return false;
  }
  return true;
}

bool Bindings::hasSites(Gene& gene, TF& tf)
{
  if (sites.find(&gene) == sites.end())
    return false;
  else
  {
    if (sites[&gene].find(&tf) == sites[&gene].end())
      return false;
  }
  return true;
}*/

void Bindings::addSite(Gene* g, TF* t, site_ptr b) 
{ 
  gene_sites_map& gsites = *(sites[g]);
  gsites[t].push_back(b); 
}

void Bindings::createScores()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      createScores(gene, tf);
    }
  }
}

void Bindings::createScores(Gene& gene)
{
  int ntfs   = tfs->size();
  
  for (int j=0; j<ntfs; j++)
  {
    TF& tf = tfs->getTF(j);
    createScores(gene, tf);
  }
}

void Bindings::createSites(ptree& pt)
{
  ptree& genes_node = pt.get_child("Genes");
  foreach_(ptree::value_type const& source, (ptree&) genes_node)
    {
      if(source.first != "Source") continue;

      foreach_(ptree::value_type& gene_node, (ptree&) source.second)
	{
	  if(gene_node.first != "Gene") continue;
	  Gene& gene = genes->getGene(gene_node.second.get<string>("<xmlattr>.name"));

	  foreach_(ptree::value_type& bindingsite_node, (ptree&) gene_node.second)
	    {
	      if(bindingsite_node.first != "BindingSite") continue;
	      TF& tf = tfs->getTF(bindingsite_node.second.get<string>("<xmlattr>.name"));

	      gene_sites_map& gsites = (*sites[&gene]);
	      site_ptr_vector& tmp_sites = gsites[&tf];
	      
	      int m = bindingsite_node.second.get<int>("<xmlattr>.m");
	      int n = bindingsite_node.second.get<int>("<xmlattr>.n");
	      char orientation = bindingsite_node.second.get<char>("<xmlattr>.orientation");
	      int nmodes = tf.getNumModes();
	      double kmax = tf.getKmax();
	      double kns = tf.getKns();
	      vector<double>& v = conc[&tf];
	      double K_exp = 0;
	      double score = 0;

	      if(mode->getBindingSiteList() == 1)
		{
		  K_exp = bindingsite_node.second.get<double>("<xmlattr>.K");
		}
	      else if(mode->getBindingSiteList() == 2)
		{
		  score = bindingsite_node.second.get<double>("<xmlattr>.score");
		  double mscore = bindingsite_node.second.get<double>("<xmlattr>.maxscore");
		  tf.setMaxScore(mscore);
		  double ddg = mscore - score;
		  double lambda = tf.getLambda();
		  
		  K_exp = exp(-ddg/lambda);
		}
	      else error("BindingSIteList must be either K or score");

	      createSite(tmp_sites, gene, tf, m, n, score, K_exp, orientation, kmax, v, nmodes, kns);
	      if(!mode->getSelfCompetition())
		trimOverlaps(gene, tf);
	      
	      updateK(gene, tf);

	      // cout << tf.getName() << " " << tmp_sites.size() << endl;
	    }
	  order_sites(gene);
	}
    }
  // printSites(cerr);
  if(mode->getVerbose() >= 2)
    cerr << "Initialized binding sites" << endl;
}

void Bindings::createSite(site_ptr_vector& tmp_sites, Gene& gene, TF& tf, int m, int n, double score, double k, char orientation, double kmax, vector<double>& v, int nmodes, double kns)
{
  site_ptr b(new BindingSite);
  b->tf = &tf;
  b->orientation = orientation;
  b->m = m;
  b->n = n;
  b->score = score;
  b->K_exp_part = k;
  double K_exp_part_times_kmax = kmax * b->K_exp_part;
  b->K_exp_part_times_kmax = K_exp_part_times_kmax;

  vector<double>& kv = b->kv;
  kv.resize(nnuc);
  for(int i = 0; i < nnuc; i++)
    {
      double tf_kmax = kmax * v[i];
      kv[i] = (b->K_exp_part * tf_kmax) / (1 + kns * tf_kmax);
    }

  b->total_occupancy.resize(nnuc);
  vector< vector<double> >& mode_occupancy = b->mode_occupancy;
  vector< vector<double> >& effective_occupancy = b->effective_occupancy;

  mode_occupancy.resize(nmodes);
  effective_occupancy.resize(nmodes);
  for (int i=0; i<nmodes; i++)
    {
      mode_occupancy[i].resize(nnuc);
      effective_occupancy[i].resize(nnuc);
    }

  b->index_in_site_map = tmp_sites.size();
  tmp_sites.push_back(b);
}


void Bindings::createSites()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    createSites(gene);
  }
}

void Bindings::createSites(Gene& gene)
{
  int ntfs   = tfs->size();

  for (int j=0; j<ntfs; j++)
  {
    TF& tf = tfs->getTF(j);
    createSites(gene, tf);
  }
  order_sites(gene);
  //for (int j=0; j<ntfs; j++)
  //{
  //  TF& tf = tfs->getTF(j);
  //  //cerr << "verify after creating 1 for " << tf.getName() << endl;
  //  // verify_order(gene, tf);
  //}
}

void Bindings::create()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    
    gene_scores_map& gscores       = *(scores[&gene]);
    gene_scores_map& saved_gscores = *(saved_scores[&gene]);
    gene_sites_map&  gsites        = *(sites[&gene]);
    gene_sites_map&  saved_gsites  = *(saved_sites[&gene]);
    
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      createScores(gene, tf);
      createSites(gene, tf);
      
      saved_gscores[&tf] = gscores[&tf];
      saved_gsites[&tf]  = gsites[&tf];
      
    }
    order_sites(gene);
    //for (int j=0; j<ntfs; j++)
    //{
    //  TF& tf = tfs->getTF(j);
    //  //cerr << "verify after creating 2 for " << tf.getName() << endl;
    //  // verify_order(gene, tf);
    //}
  }
}

void Bindings::create(Gene& gene)
{
  int ntfs   = tfs->size();

  for (int j=0; j<ntfs; j++)
  {
    TF& tf = tfs->getTF(j);
    createScores(gene, tf);
    createSites(gene, tf);
  }
  order_sites(gene);
  //for (int j=0; j<ntfs; j++)
  //  {
  //    TF& tf = tfs->getTF(j);
  //    //cerr << "verify after creating 3 for " << tf.getName() << endl;
  //    verify_order(gene, tf);
  //  }
}


// find sites for a particular gene and tf
void Bindings::createScores(Gene& gene, TF& tf) 
{
  gene_scores_map& gscores = *(scores[&gene]);
  gscores[&tf] = tf.score(gene.getSequence());
}

TFscore& Bindings::getScores(Gene& gene, TF& tf) 
{ 
  gene_scores_map& gscores = *(scores[&gene]);
  return gscores[&tf]; 
}


void Bindings::createSites(Gene& gene, TF& tf)
{
  gene_sites_map& gsites = *(sites[&gene]);
  site_ptr_vector& tmp_sites = gsites[&tf];
  
  //cerr << gene.getName() << " " << tmp_sites.size() << endl;
  if (tmp_sites.size() != 0) error("you didnt clear tmp_sites!");
  
  gene_scores_map& gscores = *(scores[&gene]);
  TFscore& t         = gscores[&tf];
  int      len       = t.mscore.size();
  int      nmodes    = tf.getNumModes();
  double   threshold = tf.getThreshold();
  double   bsize     = tf.getBindingSize();
  double   kmax      = tf.getKmax();
  double   maxscore  = tf.getMaxScore();
  double   lambda    = tf.getLambda();
  double   kns       = tf.getKns(); 
  
  //kns = kns + exp((offset-maxscore)/lambda);
  
  vector<double>& v = conc[&tf];
  
  for (int k=0; k<len; k++)
  {
    
    if (t.mscore[k] < threshold) continue;

    double fscore = t.fscore[k];
    double rscore = t.rscore[k];
    
    if (fscore >= threshold)
      createSite(tmp_sites, gene, tf, k, bsize, fscore, 'F', lambda, kmax,maxscore,v,nmodes,kns);
    if (rscore >= threshold)
      createSite(tmp_sites, gene, tf, k, bsize, rscore, 'R', lambda, kmax,maxscore,v,nmodes,kns);

  }
  if (!mode->getSelfCompetition())
    trimOverlaps(gene,tf);
  
  updateK(gene, tf);
  //printSites(tf, cerr);
}

void Bindings::add_to_ordered(Gene& gene, TF& tf)
{
  gene_sites_map& gsites       = *(sites[&gene]);
  site_ptr_vector& tf_sites    = gsites[&tf];
  vector<BindingSite*>& fsites = ordered_sites_f[&gene];
  
  int ntfsites = tf_sites.size();
  int nfsites  = fsites.size();
  int total    = nfsites + ntfsites;
  fsites.resize(total);
  
  int i, j;
  for (i=nfsites, j=0; j<ntfsites; i++, j++)
    fsites[i] = tf_sites[j].get();
}
    
  
  
  
void Bindings::createSite(site_ptr_vector& tmp_sites, Gene& gene, TF& tf,
                          int pos, double bsize, double score, char orientation, 
                          double lambda, double kmax, double maxscore,vector<double>& v, 
                          int nmodes, double kns)
{
  site_ptr b(new BindingSite);
  b->tf                    = &tf;
  b->orientation           = orientation;
  b->m                     = floor(pos  - bsize/2);
  b->n                     = b->m + bsize-1;
  b->score                 = score;
  b->pos                   = pos;
  double ddg = maxscore - b->score;

  double multiplier = 1;
  if (mode->getChromatin())
  {
    //int glength = gene.length();
    double kacc = chromatin->getKacc();
    vector<double>& acc = chromatin->getAcc(gene);
    multiplier = exp( -kacc*(1 - acc[pos]));
  }
    
  b->K_exp_part            = exp(-ddg/lambda)*multiplier;
  double K_exp_part_times_kmax = kmax * b->K_exp_part;
  b->K_exp_part_times_kmax = K_exp_part_times_kmax;
  
  vector<double>& kv = b->kv;
  kv.resize(nnuc);
  for (int i=0; i<nnuc; i++)
  {
    double tf_kmax = kmax *v[i];
    kv[i] = (b->K_exp_part * tf_kmax) / (1 + kns * tf_kmax);
  }
  
  b->total_occupancy.resize(nnuc);

  vector< vector<double> >& mode_occupancy      = b->mode_occupancy;
  vector< vector<double> >& effective_occupancy = b->effective_occupancy;
  
  mode_occupancy.resize(nmodes);
  effective_occupancy.resize(nmodes);
  for (int i=0; i<nmodes; i++)
  {
    mode_occupancy[i].resize(nnuc);
    effective_occupancy[i].resize(nnuc);
  }
  
  b->index_in_site_map = tmp_sites.size();
  tmp_sites.push_back(b);
  
  //ordered_sites_f[&gene].push_back(b.get());
}


void Bindings::trimOverlaps(Gene& gene, TF& tf)
{
  gene_sites_map& gsites = *(sites[&gene]);
  site_ptr_vector& tf_sites = gsites[&tf];
  
  int nsites = tf_sites.size();
  int i = 0;
  
  site_ptr_vector::iterator sites_begin = tf_sites.begin();
  while (i<(nsites-1))
  {
    BindingSite& site1 = *tf_sites[i];
    BindingSite& site2 = *tf_sites[i+1];
    
    if (bad_overlap_function(site1, site2))
    {
      if (site1.K_exp_part > site2.K_exp_part)
      {
        tf_sites.erase(sites_begin + i + 1);
        nsites--;
      } 
      else
      {
        tf_sites.erase(sites_begin + i);
        nsites--;
      }
    } 
    else
      i++;
  }
  
  // now we need to set the index in bindings back
  nsites = tf_sites.size();
  for (int i=0; i<nsites; i++)
  {
    BindingSite& b = *tf_sites[i];
    b.index_in_site_map = i;
  }
}
    
  
  
  
void Bindings::order_sites(Gene& gene)
{
  vector<BindingSite*>& fsites = ordered_sites_f[&gene];
  vector<BindingSite*>& rsites = ordered_sites_r[&gene];
  
  fsites.clear();
  int ntfs = tfs->size();
  gene_sites_map& gene_sites = *(sites[&gene]);
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tf_sites = gene_sites[&tf];
    int nsites = tf_sites.size();
    int cur_total = fsites.size();
    fsites.resize(cur_total+nsites);
    for (int j=0; j<nsites; j++)
    {
      fsites[cur_total + j] = tf_sites[j].get();
    }
  }
  int nsites=fsites.size();
  
  std::sort(fsites.begin(), fsites.end(), compareBindingSiteRight);

  // sort reverse sites
  nsites = fsites.size();
  
  rsites.resize(nsites);
  for (int i=0; i<nsites; i++)
    rsites[i] = fsites[nsites - i - 1];
  
  std::sort(rsites.begin(), rsites.end(), compareBindingSiteLeft);
 
  // in general, we want to know where sites are in each structure, so we will keep an index of each
  for (int i=0; i<nsites; i++)
  {
    fsites[i]->index_in_ordered_f = i;
    rsites[i]->index_in_ordered_r = i;
  }
}

void Bindings::updateK(Gene& gene, TF& tf)
{
  vector<double>&  v      = conc[&tf];
  double           kmax   = tf.getKmax();
  //double           offset = tf.getPWMOffset();
  double           kns    = tf.getKns(); 
  //double        maxscore  = tf.getMaxScore();
  //double        lambda    = tf.getLambda();
  
  //kns = kns + exp((offset-maxscore)/lambda);
  
  gene_sites_map& gsites = *(sites[&gene]);
  site_ptr_vector& tmp_sites = gsites[&tf];
  int nsites = tmp_sites.size();
  for (int k=0; k<nsites; k++)
  {
    BindingSite* b = tmp_sites[k].get();
    double K_exp_part_times_kmax = kmax * b->K_exp_part;
    b->K_exp_part_times_kmax = K_exp_part_times_kmax;
    vector<double>& kv = b->kv;
    for (int i=0; i<nnuc; i++)
    {
      double tf_kmax = kmax *v[i];
      kv[i] = (b->K_exp_part*tf_kmax) / (1 + kns*tf_kmax);
    }
  }
}

void Bindings::updateK(TF& tf)
{
  vector<double>&  v      = conc[&tf];
  double           kmax   = tf.getKmax();
  //double           offset = tf.getPWMOffset();
  double           kns    = tf.getKns(); 
  //double   maxscore  = tf.getMaxScore();
  //double   lambda    = tf.getLambda();
  
  //kns = kns + exp((offset-maxscore)/lambda);
  
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    gene_sites_map& gsites = *(sites[&gene]);
    site_ptr_vector& tmp_sites = gsites[&tf];
    int nsites = tmp_sites.size();
    for (int k=0; k<nsites; k++)
    {
      BindingSite* b = tmp_sites[k].get();
      double K_exp_part_times_kmax = kmax * b->K_exp_part;
      b->K_exp_part_times_kmax = K_exp_part_times_kmax;
      vector<double>& kv = b->kv;
      for (int i=0; i<nnuc; i++)
      {
        double tf_kmax = kmax *v[i];
        kv[i] = (b->K_exp_part*tf_kmax) / (1 + kns*tf_kmax);
      }
    }
  }
}

void Bindings::updateKandLambda(Gene& gene)
{
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    updateKandLambda(gene, tf);
  }
}

  
void Bindings::updateKandLambda(Gene& gene, TF& tf)
{
  vector<double>&     v = conc[&tf];
  
  double   kmax     = tf.getKmax();
  double   lambda   = tf.getLambda();
  double   maxscore = tf.getMaxScore();
  //double   offset   = tf.getPWMOffset();
  double   kns      = tf.getKns(); 
  
  //kns = kns + exp((offset-maxscore)/lambda);
  
  gene_sites_map& gsites = *(sites[&gene]);
  site_ptr_vector& tmp_sites = gsites[&tf];
  int nsites = tmp_sites.size();
  for (int k=0; k<nsites; k++)
  {
    BindingSite* b = tmp_sites[k].get();
    
    double ddg = maxscore - b->score;
    if (mode->getChromatin())
    {
      //int glength = gene.length();
      double kacc = chromatin->getKacc();
      vector<double>& acc = chromatin->getAcc(gene);
      ddg -= kacc * (1 - acc[b->pos]);
    }
    b->K_exp_part            = exp(-ddg/lambda);
    
    double K_exp_part_times_kmax = kmax * b->K_exp_part;
    
    b->K_exp_part_times_kmax = K_exp_part_times_kmax;
    vector<double>& kv = b->kv;
    for (int i=0; i<nnuc; i++)
    {
      double tf_kmax = kmax *v[i];
      kv[i] = (b->K_exp_part*tf_kmax) / (1 + kns*tf_kmax);
    }
  }
}

void Bindings::updateKandLambda(TF& tf)
{
  vector<double>&     v = conc[&tf];
  
  double   kmax     = tf.getKmax();
  double   lambda   = tf.getLambda();
  double   maxscore = tf.getMaxScore();
  //double   offset   = tf.getPWMOffset();
  double   kns      = tf.getKns(); 
  
  //kns = kns + exp((offset-maxscore)/lambda);
  
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    gene_sites_map& gsites = *(sites[&gene]);
    site_ptr_vector& tmp_sites = gsites[&tf];
    int nsites = tmp_sites.size();
    
    for (int k=0; k<nsites; k++)
    {
      BindingSite* b = tmp_sites[k].get();
      double ddg = maxscore - b->score;
      double multiplier = 1;
      if (mode->getChromatin())
      {
        //int glength = gene.length();
        double kacc = chromatin->getKacc();
        vector<double>& acc = chromatin->getAcc(gene);
        //cerr << acc[b->pos] << endl;
        multiplier = exp( -kacc*(1 - acc[b->pos]));
      }
      b->K_exp_part = exp(-ddg/lambda)*multiplier;
      
      double K_exp_part_times_kmax = kmax * b->K_exp_part;
      
      b->K_exp_part_times_kmax = K_exp_part_times_kmax;
      vector<double>& kv = b->kv;
      for (int i=0; i<nnuc; i++)
      {
        double tf_kmax = kmax *v[i];
        kv[i] = (b->K_exp_part*tf_kmax) / (1 + kns*tf_kmax);
      }
    }
  }
}

void Bindings::addNuc(string& nuc_id)
{
  
  //id_2_idx[nuc_id] = nnuc; 
  //idx_2_id[nnuc]   = nuc_id;
  
  IDs.push_back(nuc_id);
  nnuc++;
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    const string& tfname = tf.getName();
    conc[&tf].push_back(tfdata->getDataPoint("TF",tfname, "ID",nuc_id));
  }
}



void Bindings::saveScores(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    saveScores(gene, tf);
  }
}

void Bindings::saveScores(Gene& gene, TF& tf)
{
  gene_scores_map& gscores       = *(scores[&gene]);
  gene_scores_map& saved_gscores = *(saved_scores[&gene]);

  saved_gscores[&tf] = gscores[&tf];

}

void Bindings::restoreScores(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    restoreScores(gene, tf);
  }
}

void Bindings::restoreScores(Gene& gene, TF& tf)
{
  gene_scores_map& gscores       = *(scores[&gene]);
  gene_scores_map& saved_gscores = *(saved_scores[&gene]);
  gscores[&tf] = saved_gscores[&tf];
}

void Bindings::updateScores(Gene& gene, TF& tf)
{
  gene_scores_map& gscores = *(scores[&gene]);
  gscores[&tf] = tf.score(gene.getSequence());
}

void Bindings::updateScores(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    updateScores(gene, tf);
  }
}

void Bindings::updateScores()
{
  int ntfs   = tfs->size();
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs->getTF(j);
      updateScores(gene, tf);
    }
  }
}

void Bindings::updateScores(Gene& gene)
{
  gene_scores_map& gscores = *(scores[&gene]);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    gscores[&tf] = tf.score(gene.getSequence());
  }
}

void Bindings::updateSites()
{
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    updateSites(tf);
  }
}

void Bindings::updateSites(Gene& gene)
{
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    updateSites(gene, tf);
  }
}

void Bindings::updateSites(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    updateSites(gene,tf);
  }
}

void Bindings::eraseTF(Gene& gene, TF& tf)
{
  /*
  vector<BindingSite*>& fsites = ordered_sites_f[&gene];
  vector<BindingSite*>& rsites = ordered_sites_r[&gene];
  
  // these are populated by base, so they are in order!
  site_ptr_vector& tf_sites = sites[&gene][&tf];
  
  vector<BindingSite*>::iterator fbegin = fsites.begin();
  vector<BindingSite*>::iterator rbegin = rsites.begin();
  
  int ntfsites = tf_sites.size();
  int i, j;
  for (i=0, j=(ntfsites-1); i<ntfsites; i++, j--)
  {
    BindingSite& fsite = *tf_sites[j];
    BindingSite& rsite = *tf_sites[i];

    fsites.erase(fbegin + fsite.index_in_ordered_f);
    rsites.erase(rbegin + rsite.index_in_ordered_r);
  }
  */
  gene_sites_map& gsites = *(sites[&gene]);
  site_ptr_vector& tf_sites = gsites[&tf];
  tf_sites.clear();
} 
  

  
void Bindings::updateSites(Gene& gene, TF& tf)
{
  //cerr << "verify before update for " << tf.getName() << endl;
  //verify_order(gene, tf);
  eraseTF(gene,tf);
  createSites(gene,tf);
  order_sites(gene);
  //cerr << "verify after update for " << tf.getName() << endl;
  //verify_order(gene, tf);
}

/* 
a function to verify that the ordered sites and tfs all point to the right 
places. We need to do do comparisons to do wo. First, we need to verify that 
everything in the site_ptr_vector is in the ordered sites at the right place,
second, we need to verify that all the ordered sites point to valid objects that
are in the site_ptr_vector
*/

void Bindings::verify_order(Gene& gene, TF& tf)
{
  vector<BindingSite*>& fsites = ordered_sites_f[&gene];                       
  vector<BindingSite*>& rsites = ordered_sites_r[&gene];
  gene_sites_map& gsites = *(sites[&gene]);
  site_ptr_vector& tf_sites    = gsites[&tf];
  
  // verify that ever site points to something in the ordered list
  int nsites = tf_sites.size();
  //cerr << nsites << " for tf " << tf.getName() << endl;
  for (int j=0; j<nsites; j++)
  {
    BindingSite& b = *tf_sites[j];
    // check that this index exists in the sites
    if (b.index_in_ordered_f >= (int) fsites.size())
      error("Tried to access index " + to_string_(b.index_in_ordered_f) + " of " + to_string_(fsites.size()) + " from forward ordered sites");
    if (b.index_in_ordered_r >= (int) rsites.size())
      error("Tried to access index " + to_string_(b.index_in_ordered_r) + " of " + to_string_(rsites.size()) + " from reverse ordered sites");
    // check that this points to itself
    if (b.index_in_site_map != j)
      error(" index in site map malformed for site " + to_string_(j) + " of tf " + tf.getName());
    
    if (fsites[b.index_in_ordered_f] != &b)
      error(" index in ordered forward malformed for site " + to_string_(j) + " of tf " + tf.getName());
    
    if (rsites[b.index_in_ordered_r] != &b)
      error(" index in ordered reverse malformed for site " + to_string_(j) + " of tf " + tf.getName());
  }
  
  // verify that the ordered sites are formed correctly
  nsites = fsites.size();
  if (nsites != (int) rsites.size())
    error("different numbers of sites in forward and reverse ordered lists");
  
  for (int i=0; i<nsites; i++)
  {
    BindingSite& sitef = *fsites[i];
    BindingSite& siter = *rsites[i];
    
    TF& tff = *sitef.tf;
    TF& tfr = *siter.tf;
    
    if (&tff == &tf)
    {
      if (sitef.index_in_site_map >= (int) tf_sites.size())
        error("Tried to access index " + to_string_(sitef.index_in_site_map) + " of " + to_string_(tf_sites.size()) + " from sites for tf " + tff.getName());
      if (tf_sites[sitef.index_in_site_map].get() != &sitef)
        error(" ordered f site " + to_string_(i) + " does not point to a valid site" + to_string_(sitef.index_in_ordered_f));
    }
    if (&tfr == &tf)
    {
      if (siter.index_in_site_map >= (int) tf_sites.size())
        error("Tried to access index " + to_string_(siter.index_in_site_map) + " of " + to_string_(tf_sites.size()) + " from sites for tf " + tfr.getName());
      if (tf_sites[siter.index_in_site_map].get() != &siter)
        error(" ordered r site " + to_string_(i) + " does not point to a valid site " + to_string_(siter.index_in_ordered_r));
    }
  }
}
    
    

void::Bindings::saveSites(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    saveSites(gene, tf);
  }
}

void::Bindings::saveSites(Gene& gene, TF& tf)
{  
  gene_sites_map& saved_gsites = *(saved_sites[&gene]);
  gene_sites_map& gsites       = *(sites[&gene]);
  
  saved_gsites[&tf] = gsites[&tf];
}

void::Bindings::restoreSites(TF& tf)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    restoreSites(gene, tf);
  }
    
}

void::Bindings::restoreSites(Gene& gene, TF& tf)
{
  //for (int i=0; i<tfs->size(); i++)
  //{
  //  TF& tf_test = tfs->getTF(i);
  // // printSites(gene,tf_test,cerr);
  //  //cerr << "verify before restore for " << tf_test.getName() << endl;
  //  verify_order(gene, tf_test);
  //}
  gene_sites_map& saved_gsites = *(saved_sites[&gene]);
  gene_sites_map& gsites       = *(sites[&gene]);
  gsites[&tf] = saved_gsites[&tf];
  //sites[&gene][&tf] = saved_sites[&gene][&tf];
  order_sites(gene); 
  //cerr << "fsites.size() = " << ordered_sites_f[&gene].size() << endl;
  //cerr << "rsites.size() = " << ordered_sites_r[&gene].size() << endl;
  
  //for (int i=0; i<tfs->size(); i++)
  //{
  //  TF& tf_test = tfs->getTF(i);
  //  //printSites(gene,tf_test,cerr);
  //  //cerr << "verify after restore for " << tf_test.getName() << endl;
  //  verify_order(gene, tf_test);
  //}
}


void Bindings::saveOccupancy()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
    saveOccupancy(genes->getGene(i));
}


void Bindings::saveOccupancy(Gene& gene)
{
  gene_sites_map& gsites = *(sites[&gene]);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
      tfsites[j]->saved_total_occupancy = tfsites[j]->total_occupancy;
  }
}
        

void Bindings::restoreOccupancy()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
    restoreOccupancy(genes->getGene(i));
}

void Bindings::restoreOccupancy(Gene& gene)
{
  gene_sites_map& gsites = *(sites[&gene]);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
      tfsites[j]->total_occupancy = tfsites[j]->saved_total_occupancy;
  }
}


/* this function saves effective occupancy and resets occupancy to 
simple fractional occupancy */
void Bindings::saveEffectiveOccupancy()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
    saveEffectiveOccupancy(genes->getGene(i));
}

void Bindings::saveEffectiveOccupancy(Gene& gene)
{
  gene_sites_map& gsites = *(sites[&gene]);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
      tfsites[j]->saved_effective_occupancy = tfsites[j]->effective_occupancy;
  }
}

void Bindings::restoreEffectiveOccupancy()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    restoreEffectiveOccupancy(gene);
  }
}

void Bindings::restoreEffectiveOccupancy(Gene& gene)
{
  gene_sites_map& gsites = *(sites[&gene]);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
      tfsites[j]->effective_occupancy = tfsites[j]->saved_effective_occupancy;
  }
}

void Bindings::saveModeOccupancy()
{
  int ngenes = genes->size();
  
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    saveModeOccupancy(gene);
  }
}

void Bindings::saveModeOccupancy(Gene& gene)
{
  gene_sites_map& gsites = *(sites[&gene]);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
    {
      tfsites[j]->saved_effective_occupancy = tfsites[j]->effective_occupancy;
      tfsites[j]->saved_mode_occupancy      = tfsites[j]->mode_occupancy;
    }
  }
}


void Bindings::restoreModeOccupancy()
{
  int ngenes = genes->size();

  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    restoreModeOccupancy(gene);
  }
}

void Bindings::restoreModeOccupancy(Gene& gene)
{
  gene_sites_map& gsites = *(sites[&gene]);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
    {
      tfsites[j]->effective_occupancy = tfsites[j]->saved_effective_occupancy;
      tfsites[j]->mode_occupancy      = tfsites[j]->saved_mode_occupancy;
    }
  }
}


void Bindings::saveAllOccupancy()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
    saveAllOccupancy(genes->getGene(i));
}

void Bindings::saveAllOccupancy(Gene& gene)
{
  gene_sites_map& gsites = *(sites[&gene]);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
    {
      site_ptr b = tfsites[j];
      b->saved_kv                  = b->kv;
      b->saved_effective_occupancy = b->effective_occupancy;
      b->saved_mode_occupancy      = b->mode_occupancy;
      b->saved_total_occupancy     = b->total_occupancy;
    }
  }
}

void Bindings::restoreAllOccupancy()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
    restoreAllOccupancy(genes->getGene(i));
}

void Bindings::restoreAllOccupancy(Gene& gene)
{
  gene_sites_map& gsites = *(sites[&gene]);
  
  int ntfs = tfs->size();
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    site_ptr_vector& tfsites = gsites[&tf];
    int nsites = tfsites.size();
    
    for (int j=0; j<nsites; j++)
    {
      site_ptr b = tfsites[j];
      b->kv                  = b->saved_kv;
      b->mode_occupancy      = b->saved_mode_occupancy;
      b->effective_occupancy = b->saved_effective_occupancy;
      b->total_occupancy     = b->saved_total_occupancy;
    }
  }
}


bool Bindings::isEqual(Bindings& test_bindings)
{
  // make sure we have the same number of tfs and genes
  genes_ptr test_genes = test_bindings.getGenes();
  tfs_ptr   test_tfs   = test_bindings.getTFs();
  
  int n_test_genes = test_genes->size();
  int n_test_tfs   = test_tfs->size();
  
  int n_genes = genes->size();
  int n_tfs   = tfs->size();
  
  if (n_test_genes != n_genes) warning("Number of genes not equal!");
  if (n_test_tfs != n_tfs) warning("Number of tfs not equal!");
  
  // make sure they are the same tfs and genes
  for (int i=0; i<n_genes; i++)
  {
    Gene& test_gene = test_genes->getGene(i);
    Gene& gene      = genes->getGene(i);
    if (test_gene.getName() != gene.getName()) warning("Genes not the same!");
    for (int j=0; j<n_tfs; j++)
    {
      TF& test_tf = test_tfs->getTF(j);
      TF& tf      = tfs->getTF(j);
      if (test_tf.getName() != tf.getName()) warning("TFs not the same!");
      TFscore& test_score = test_bindings.getScores(test_gene, test_tf);
      TFscore& score      = getScores(gene, tf);
      int n = test_score.fscore.size();
      for (int k=0; k<n; k++)
      {
        if (test_score.fscore[k] != score.fscore[k]) error("fscores not equal");
        if (test_score.rscore[k] != score.rscore[k]) error("rscores not equal");
      }
    }
    //vector<BindingSite*>& test_sites_f = test_bindings.getFsites(test_gene);
    //vector<BindingSite*>& test_sites_r = test_bindings.getRsites(test_gene);
    //
    //int n = ordered_sites_f.size();
    //for (int j=0; j<n; j++)
    //{
    //}
      
  }
    
  return true;
}
  

void Bindings::printSites(ostream& os)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    printSites(gene, os);
  }
}

void Bindings::printSites(Gene& gene, ostream& os)
{
  os << gene.getName() << endl;
  int print_source = 0;
  int ntfs = tfs->size();
  
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    const string& pwm_source = tf.getPWMSource();
    if (pwm_source != string(""))
      print_source = max(print_source, (int) pwm_source.size());
  }
  
  int p = mode->getPrecision();
  int w = p+7;
  os << setprecision(p)
     << setw(w) << "name"
     << setw(w) << "index"
     << setw(w) << "index_f"
     << setw(w) << "index_r"
     << setw(w) << "start"
     << setw(w) << "end"
     << setw(w) << "score"
     << setw(w) << "maxscore"
     << setw(w) << "K"
     << setw(w) << "K*kmax"
     << setw(w) << "strand";
     
  if (print_source)
    os << setw(print_source+2) << "source";
  
  os << endl;
  for (int i=0; i<ntfs; i++)
    printSites(gene, tfs->getTF(i), os, p, print_source);
}

void Bindings::printSites(TF& tf, ostream& os)
{

  int print_source = 0;
  const string& pwm_source = tf.getPWMSource();
  if (pwm_source != string(""))
    print_source = max(print_source, (int) pwm_source.size());
  
  int p = mode->getPrecision();
  int w = p+7;
  
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    os << gene.getName() << endl;
    os << setprecision(p)
     << setw(w) << "name"
     << setw(w) << "index"
     << setw(w) << "index_f"
     << setw(w) << "index_r"
     << setw(w) << "start"
     << setw(w) << "end"
     << setw(w) << "score"
     << setw(w) << "K"
     << setw(w) << "K*kmax"
     << setw(w) << "strand";
     
    if (print_source)
      os << setw(print_source+2) << "source";
    
    os << endl;
    printSites(gene, tf, os, p, print_source);
  }
}

void Bindings::printSites(Gene& gene, TF& tf, ostream& os)
{
  int print_source = 0;
  int ntfs = tfs->size();
  
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    const string& pwm_source = tf.getPWMSource();
    if (pwm_source != string(""))
      print_source = max(print_source, (int) pwm_source.size());
  }
  
  int p = mode->getPrecision();
  int w = p+7;
  os << setprecision(p)
     << setw(w) << "name"
     << setw(w) << "index"
     << setw(w) << "index_f"
     << setw(w) << "index_r"
     << setw(w) << "start"
     << setw(w) << "end"
     << setw(w) << "score"
     << setw(w) << "K"
     << setw(w) << "K*kmax"
     << setw(w) << "strand";
     
  if (print_source)
    os << setw(print_source+2) << "source";
  
  os << endl;
  printSites(gene, tf, os, p, print_source);
}
  
  
void Bindings::printSites(Gene& gene, TF& tf, ostream& os, int p, int print_source)
{ 
  int w = p+7;
  gene_sites_map& gsites = *(sites[&gene]);
  site_ptr_vector& tmp_sites = gsites[&tf];
  //double maxscore   = tf.getMaxScore();
  int    left_bound = gene.getLeftBound();
  int    nsites     = tmp_sites.size();
  for (int i=0; i<nsites; i++)
  {
    BindingSite& b = *tmp_sites[i];
    os << setprecision(p)
       << setw(w) << b.tf->getName()
       << setw(w) << b.index_in_site_map
       << setw(w) << b.index_in_ordered_f
       << setw(w) << b.index_in_ordered_r
       << setw(w) << b.m + left_bound
       << setw(w) << b.n + left_bound
       << setw(w) << b.score
       << setw(w) << b.tf->getMaxScore()
       << setw(w) << b.K_exp_part
       << setw(w) << b.K_exp_part_times_kmax;
     if (b.orientation == 'F')
       os << setw(w) << "+";
     else 
       os << setw(w) << "-";
     
     if (print_source)
        os << setw(print_source+2) << tmp_sites[i]->tf->getPWMSource();
      os << endl;
  }
}

void Bindings::printScores(Gene& gene, ostream& os)
{
  int p = mode->getPrecision();
  int w = p+7;
  
  gene_scores_map& gscores = *(scores[&gene]);
  
  int ntfs = tfs->size();
  
  os << setw(w) << "base";
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    os << setw(w) << tf.getName() + "_f";
    os << setw(w) << tf.getName() + "_r";
  }
  os << endl;
  
  int length     = gene.length();
  int left_bound = gene.getLeftBound();
  
  for (int bp=0; bp<length; bp++)
  {
    os << setw(w) << bp + left_bound;
    for (int i=0; i<ntfs; i++)
    {
      TF& tf = tfs->getTF(i);
      os << setw(w) << gscores[&tf].fscore[bp];
      os << setw(w) << gscores[&tf].rscore[bp];
    }
    os << endl;
  }
}
      

void Bindings::printTotalOccupancy(Gene& gene, ostream& os, bool invert)
{
  int p = mode->getPrecision();
  int w = p + 7;
  os << setprecision(p);
  
  gene_sites_map& gsites = *(sites[&gene]);
  int ntfs = tfs->size();
  if (invert == false)
  {
    // loop through everything and figure out the header
    os << setw(w) << "id";
  
    // find all the ids, using map to prevent duplicates
    // simultaneously print binding site headers
    for (int i = 0; i<ntfs; i++)
    {
      TF& tf = tfs->getTF(i);
      site_ptr_vector& tmp_sites = gsites[&tf];
      int nsites = tmp_sites.size();
      for (int i=0; i<nsites; i++)
      {
        string header = tf.getName();
        stringstream convert;
        convert << i;
        string num = convert.str();
        header += num;
        //header << i;
        os << setw(w) <<  header;
      }
    }
  
    os << endl;
  
    for (int nuc=0; nuc<nnuc; nuc++)
    {
      string& id = IDs[nuc];
      os << setw(w) << id;
    
      for (int i = 0; i<ntfs; i++)
      {
        TF& tf = tfs->getTF(i);
        site_ptr_vector& tmp_sites = gsites[&tf];
        int nsites = tmp_sites.size();
        for (int i=0; i<nsites; i++)
          os << setw(w) << tmp_sites[i]->total_occupancy[nuc];
      }
      os << endl;
    }
  } 
  else 
  {
    os << setw(w) << "site";
    for (int nuc=0; nuc<nnuc; nuc++)
    {
      string& id = IDs[nuc];
      os << setw(w) << id;
    }
    os << endl;
    for (int i = 0; i<ntfs; i++)
    {
      TF& tf = tfs->getTF(i);
      site_ptr_vector& tmp_sites = gsites[&tf];
      int nsites = tmp_sites.size();
      for (int i=0; i<nsites; i++)
      {
        string header = tf.getName();
        stringstream convert;
        convert << i;
        string num = convert.str();
        header += num;
        //header << i;
        os << setw(w) <<  header;
        for (int nuc=0; nuc<nnuc; nuc++)
        {
          os << setw(w) << tmp_sites[i]->total_occupancy[nuc];
        }
      os << endl;
      }
    }
  }
}

void Bindings::printEffectiveOccupancy(Gene& gene, ostream& os, bool invert)
{
  int w = 12;
  os << setprecision(3);
  
  gene_sites_map& gsites = *(sites[&gene]);
  int ntfs = tfs->size();
  //int ngenes = genes->size();
  if (!invert)
  {
    
    // loop through everything and figure out the header
    os << setw(w) << "id";
  
    // find all the ids, using map to prevent duplicates
    // simultaneously print binding site headers
    for (int k=0; k<ntfs; k++)
    {
      TF& tf = tfs->getTF(k);
      site_ptr_vector& tmp_sites = gsites[&tf];
      int nsites = tmp_sites.size();
      int nmodes = tf.getNumModes();
      for (int i=0; i<nsites; i++)
      {      
        for (int j=0; j<nmodes; j++)
        {
          string header = tf.getName();
          stringstream convert;
          convert << i;
          if (j>0) convert << "_" << j+1;
          string num = convert.str();
          header += num;
          //header << i;
          os << setw(w) <<  header;
        }
      }
    }
  
    os << endl;
  
    for (int nuc=0; nuc<nnuc; nuc++)
    {
      string& id = IDs[nuc];
      os << setw(w) << id;
   
      for (int k=0; k<ntfs; k++)
      {
        TF& tf = tfs->getTF(k);
        site_ptr_vector& tmp_sites = gsites[&tf];
        int nsites = tmp_sites.size();
        int nmodes = tf.getNumModes();
        for (int i=0; i<nsites; i++)
        {
          for (int j=0; j<nmodes; j++)
            os << setw(w) << tmp_sites[i]->effective_occupancy[j][nuc];
        }
      }
      os << endl;
    }
  }
  else 
  {
    os << setw(w) << "site";
    for (int nuc=0; nuc<nnuc; nuc++)
    {
      string& id = IDs[nuc];
      os << setw(w) << id;
    }
    os << endl;

    for (int k=0; k<ntfs; k++)
    {
      TF& tf = tfs->getTF(k);
      site_ptr_vector& tmp_sites = gsites[&tf];
      int nsites = tmp_sites.size();
      int nmodes = tf.getNumModes();
      for (int i=0; i<nsites; i++)
      {
        for (int j=0; j<nmodes; j++)
        {
          string header = tf.getName();
          stringstream convert;
          convert << i;
          if (j>0) convert << "_" << j+1;
          string num = convert.str();
          header += num;
          //header << i;
          os << setw(w) <<  header;
          for (int nuc=0; nuc<nnuc; nuc++)
          {
            os << setw(w) << tmp_sites[i]->effective_occupancy[j][nuc];
          }
        os << endl;
        }
      }
    }
  }
}

void Bindings::printModeOccupancy(Gene& gene, ostream& os, bool invert)
{
  int w = 12;
  os << setprecision(3);
  int ntfs = tfs->size();
  //int ngenes = genes->size();
  
  gene_sites_map& gsites = *(sites[&gene]);
  
  if (!invert)
  {
    
    // loop through everything and figure out the header
    os << setw(w) << "id";
  
    // find all the ids, using map to prevent duplicates
    // simultaneously print binding site headers
    for (int k=0; k<ntfs; k++)
    {
      TF& tf = tfs->getTF(k);
      site_ptr_vector& tmp_sites = gsites[&tf];
      int nsites = tmp_sites.size();
      int nmodes = tf.getNumModes();
      for (int i=0; i<nsites; i++)
      {      
        for (int j=0; j<nmodes; j++)
        {
          string header = tf.getName();
          stringstream convert;
          convert << i;
          if (j>0) convert << "_" << j+1;
          string num = convert.str();
          header += num;
          //header << i;
          os << setw(w) <<  header;
        }
      }
    }
  
    os << endl;
  
    for (int nuc=0; nuc<nnuc; nuc++)
    {
      string& id = IDs[nuc];
      os << setw(w) << id;
   
      for (int k=0; k<ntfs; k++)
      {
        TF& tf = tfs->getTF(k);
        site_ptr_vector& tmp_sites = gsites[&tf];
        int nsites = tmp_sites.size();
        int nmodes = tf.getNumModes();
        for (int i=0; i<nsites; i++)
        {
          for (int j=0; j<nmodes; j++)
            os << setw(w) << tmp_sites[i]->mode_occupancy[j][nuc];
        }
      }
      os << endl;
    }
  }
  else 
  {
    os << setw(w) << "site";
    for (int nuc=0; nuc<nnuc; nuc++)
    {
      string& id = IDs[nuc];
      os << setw(w) << id;
    }
    os << endl;

    for (int k=0; k<ntfs; k++)
    {
      TF& tf = tfs->getTF(k);
      site_ptr_vector& tmp_sites = gsites[&tf];
      int nsites = tmp_sites.size();
      int nmodes = tf.getNumModes();
      for (int i=0; i<nsites; i++)
      {
        for (int j=0; j<nmodes; j++)
        {
          string header = tf.getName();
          stringstream convert;
          convert << i;
          if (j>0) convert << "_" << j+1;
          string num = convert.str();
          header += num;
          //header << i;
          os << setw(w) <<  header;
          for (int nuc=0; nuc<nnuc; nuc++)
          {
            os << setw(w) << tmp_sites[i]->mode_occupancy[j][nuc];
          }
        os << endl;
        }
      }
    }
  }
}

void Bindings::write(ptree& output)
{
  ptree& genes_out = output.get_child("Genes");
  
  foreach_(ptree::value_type& source_out, genes_out)
    {
      if(source_out.first != "Source") continue;

      foreach_(ptree::value_type& gene_out, (ptree&) source_out.second)
	{
	  if(gene_out.first != "Gene") continue;

	  int ntfs = tfs->size();
	  Gene& gene = genes->getGene(gene_out.second.get<string>("<xmlattr>.name"));
	  if(!gene.getInclude()) continue;
	  
	  for(int i = 0; i < ntfs; i++)
	    {
	      TF& tf = tfs->getTF(i);
	      gene_sites_map& gsites = (*sites[&gene]);
	      site_ptr_vector& tmp_sites = gsites[&tf];
	      int nsites = tmp_sites.size();

	      // cout << "TF " << tf.getName() << "\t" << nsites << endl;
	      for(int j = 0; j < nsites; j++)
		{
		  write((ptree&) gene_out.second, *tmp_sites[j]);
		}
	    }
	}
    }
}

void Bindings::write(ptree& pt, BindingSite& b) const
{
  ptree& bsite = pt.add("BindingSite", "");
  bsite.put("<xmlattr>.name", b.tf->getName());
  bsite.put("<xmlattr>.m", b.m);
  bsite.put("<xmlattr>.n", b.n);
  bsite.put("<xmlattr>.orientation", b.orientation);
  bsite.put("<xmlattr>.score", b.score);
  bsite.put("<xmlattr>.maxscore", b.tf->getMaxScore());
  bsite.put("<xmlattr>.K", b.K_exp_part);
}
