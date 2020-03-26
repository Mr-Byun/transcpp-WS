/*********************************************************************************
*                                                                                *
*     quenching.cpp                                                              *
*                                                                                *
*     Contains structure and methods to hold quenching interactions              *
*                                                                                *
*********************************************************************************/

#include "quenching.h"
#include <boost/foreach.hpp>
#include <limits>

# define foreach_ BOOST_FOREACH

QuenchingInteractions::QuenchingInteractions() {}

void QuenchingInteractions
::create(genes_ptr g, tfs_ptr t, bindings_ptr b, distances_ptr d)
{
  genes     = g;
  tfs       = t;
  bindings  = b;
  distances = d;
  
  dist = distances->getDistance("Quenching");
  
  int ngenes = genes->size();
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    
    gene_quenches_ptr gquenches(new gene_quenches);
    gene_quenches_ptr saved_gquenches(new gene_quenches);
    
    quenches[&gene]       = gquenches;
    saved_quenches[&gene] = saved_gquenches;
    
    int ntfs   = tfs->size();
    for (int i = 0; i<ntfs; i++) // loop through actors
    {
      TF& tf1 = tfs->getTF(i);
      if (tf1.neverQuenches()) continue; // skip if never a quencher
      for (int j=0; j<ntfs; j++) // loop through targets
      {
        TF& tf2 = tfs->getTF(j);
        if (tf2.neverActivates()) continue; // skip if never an activator
        set(gene, tf1, tf2);
      }
    }   
  }
  //printSummary();
}

void QuenchingInteractions
::set(Gene& gene, TF& actor, TF& target)
{
  site_ptr_vector::iterator aiter;
  site_ptr_vector::iterator titer;
  
  site_ptr_vector& actorsites  = bindings->getSites(gene, actor);
  site_ptr_vector& targetsites = bindings->getSites(gene, target);
  
  gene_quenches& gquenches = *(quenches[&gene]);
  vector<QuenchingInteraction>& q = gquenches[&actor][&target];

  double max_dist     = dist->getMaxDistance();
  //int    nactors      = actorsites.size();
  int    ntargets     = targetsites.size();
  int    start = 0;
  
  for (aiter=actorsites.begin(); aiter != actorsites.end(); ++aiter)
  {
    int c = 0; // value used to see if we 
    BindingSite& actor     = *(*aiter);
    //site_ptr     actor_ptr = *aiter;
    int m1 = actor.m;
    int n1 = actor.n;
    
    // note this only works if sites are ordered! (which they are because we did a linear search on DNA)
    for (int i=start; (i<ntargets && c!=2); i++)
    {
      site_ptr     target_ptr = targetsites[i];
      BindingSite& target     = *target_ptr;
      
      int m2 = target.m;
      int n2 = target.n;
      double d;
      
      int dm = abs(m1 - n2);
      int dn = abs(n1 - m2);
      
      if (dn <= dm)
        d = dn;
      else
        d = dm;
      
      bool  overlapped = (m1 < n2 && m2 < n1);
      double df        = dist->getDistFunc(d);
      if ( d < max_dist)
      {
        if (c==0) // we found the first site in range
        {
          c++; 
          start = i; // since this is the first site in range, and the next site is farther down, this is the starting point for our next loop
        }
        if (!overlapped && df > 0)
        {
          QuenchingInteraction quench;
          quench.actor  = *aiter;
          quench.target = target_ptr;
          quench.distcoef = df;
          
          q.push_back(quench);
        }
      } 
      else
      {
        if (c==1) c++; // we found the last site in range
      }
    }
  }
}
  

/* I this function broke when I made this safe for parallelization over genes,
but it is never used, so I just commented it out 

bool QuenchingInteractions
::hasQuenchingInteractions(Gene& gene, TF& actor, TF& target)
{
  if (quenches.find(&gene) == quenches.end())
    return false;
  else
  {
    if (quenches[&gene].find(&actor) == quenches[&gene].end())
      return false;
    else
    {
      if (quenches[&gene][&actor].find(&target) == quenches[&gene][&actor].end())
        return false;
    }
  }
  return true;
}
*/

void QuenchingInteractions
::calc()
{
  int ngenes = genes->size();
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    calc(gene);
  }
}

void QuenchingInteractions
::calc(Gene& gene)
{
  initialize(gene);
    
  int ntfs   = tfs->size();
  for (int i = 0; i<ntfs; i++) // loop through actors
  {
    TF& tf1 = tfs->getTF(i);
    if (tf1.neverQuenches()) continue; // skip if never a quencher
    for (int j=0; j<ntfs; j++) // loop through targets
    {
      TF& tf2 = tfs->getTF(j);
      if (tf2.neverActivates()) continue; // skip if never an activator
      calc(gene, tf1, tf2);
    }
  } 
}
  
void QuenchingInteractions
::calc(Gene& gene, TF& actor, TF& target)
{
  gene_quenches& gquenches = *(quenches[&gene]);
  vector<QuenchingInteraction>& quench_vector = gquenches[&actor][&target];
  
  vector<double> actor_coefs  = actor.getCoefs();
  vector<double> target_coefs = target.getCoefs();

  int n_actor_modes  = actor_coefs.size();
  int n_target_modes = target_coefs.size();
  
  int nquench = quench_vector.size();
  
  for (int i=0; i<nquench; i++)
  {
    QuenchingInteraction& quench      = quench_vector[i];
    BindingSite& actor_site  = *quench.actor;
    BindingSite& target_site = *quench.target;
    double distcoef          = quench.distcoef;

    vector<double>* actor_occupancy;
    vector<double>* target_occupancy;
    
    for (int j=0; j<n_actor_modes; j++)
    {
      double efficiency = -actor_coefs[j];
      if (efficiency <= 0) continue;
      double efd = efficiency * distcoef;
      actor_occupancy = &(actor_site.mode_occupancy[j]);
      for (int k=0; k<n_target_modes; k++)
      {
        if (target_coefs[k] <= 0) continue;
        
        target_occupancy = &(target_site.effective_occupancy[k]);
        quench_f(*actor_occupancy, *target_occupancy, efd);
      }
    }
  }
}


void QuenchingInteractions
::quench_f(vector<double>& actor_vec, vector<double>& target_vec, double efd)
{
  vector<double>::iterator i;
  vector<double>::iterator j;
  //int n = actor_vec.size();
  for (i = actor_vec.begin(), j = target_vec.begin() ; i != actor_vec.end(); ++i, ++j)
  {
    double reduction = 1 - (*i) * efd;
    (*j) *= reduction;
  }
}

void QuenchingInteractions
::save()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    save(gene);
  }
}

void QuenchingInteractions
::save(Gene& gene)
{
  gene_quenches& gquenches       = *(quenches[&gene]);
  gene_quenches& saved_gquenches = *(saved_quenches[&gene]);
  saved_gquenches = gquenches;
}

void QuenchingInteractions
::clear()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    clear(gene);
  }
}

void QuenchingInteractions
::clear(Gene& gene)
{
  gene_quenches& gquenches = *(quenches[&gene]);
  gquenches.clear();
}

void QuenchingInteractions
::restore()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    restore(gene);
  }
}

void QuenchingInteractions
::restore(Gene& gene)
{
  gene_quenches& gquenches       = *(quenches[&gene]);
  gene_quenches& saved_gquenches = *(saved_quenches[&gene]);
  gquenches = saved_gquenches;
}

void QuenchingInteractions
::update()
{
  int ngenes = genes->size();
 
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    update(gene);
  }
}

void QuenchingInteractions
::update(Gene& gene)
{
  int ntfs   = tfs->size();

  gene_quenches& gquenches       = *(quenches[&gene]);
  
  for (int i = 0; i<ntfs; i++) // loop through actors
  {
    TF& tf1 = tfs->getTF(i);
    if (tf1.neverQuenches())
    {
      gquenches[&tf1].clear(); // make sure it is empty if it never quenches
      
      continue;                      // skip if not a quencher
    }
    for (int j=0; j<ntfs; j++) // loop through targets
    {
      TF& tf2 = tfs->getTF(j);
      gquenches[&tf1][&tf2].clear(); // make sure it is empty
      if (tf2.neverActivates()) continue;  // skip if never an activator
      set(gene, tf1, tf2);
    }
  } 
}
  

void QuenchingInteractions
::printSummary()
{
  int ngenes = genes->size();
  int ntfs   = tfs->size();
  int total  = 0;
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    gene_quenches& gquenches = *(quenches[&gene]);
    cerr << gene.getName() << endl;
    for (int j=0; j<ntfs; j++)
    {
      TF& tf1 = tfs->getTF(j);
      for (int k=0; k<ntfs; k++)
      {
        TF& tf2 = tfs->getTF(k);
        int nq = gquenches[&tf1][&tf2].size();
        cerr << tf1.getName() << " -> " << tf2.getName() << " : " << nq << endl;
        total += nq;
      }
    }
  }
  cerr << "Total: " << total << endl << endl;
}

void QuenchingInteractions
::initialize()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    initialize(gene);
  }
}

void QuenchingInteractions
::initialize(Gene& gene)
{
  int ntfs   = tfs->size();
  for (int j = 0; j<ntfs; j++) 
  {
    TF& tf = tfs->getTF(j);
    initialize(gene, tf);
  }
}
      

void QuenchingInteractions
::initialize(Gene& gene, TF& tf)
{
  site_ptr_vector& sites = bindings->getSites(gene, tf);

  vector<double> coefs = tf.getCoefs();
  
  int nmodes = coefs.size();
  int nsites = sites.size();

  for (int i=0; i<nsites; i++)
  {
    BindingSite* site = sites[i].get();
    int nnuc = site->total_occupancy.size();
    for (int j=0; j<nmodes; j++)
    {
      if (coefs[j] > 0)
        site->effective_occupancy[j] = site->mode_occupancy[j];
      else
      {
        for (int k=0; k<nnuc; k++)
          site->effective_occupancy[j][k] = 0;
      }
    }
  }
}




ModifyingInteractions::ModifyingInteractions() {}


void ModifyingInteractions
::create(genes_ptr g, tfs_ptr t, bindings_ptr b, distances_ptr d)
{
  genes     = g;
  tfs       = t;
  bindings  = b;
  distances = d;
  
  dist = distances->getDistance("Quenching");

  int ngenes = genes->size();
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    
    gene_mods_ptr gmods(new gene_mods);
    gene_mods_ptr saved_gmods(new gene_mods);
    
    mods[&gene]       = gmods;
    saved_mods[&gene] = saved_gmods;
    
    int ntfs   = tfs->size();
    for (int i = 0; i<ntfs; i++) // loop through actors
    {
      TF& tf1 = tfs->getTF(i);
      vector< pair<TF*, coeffect_ptr> >& targets = tf1.getTargets();
      int ntargets = targets.size();
      for (int j=0; j<ntargets; j++)
      {
        pair<TF*, coeffect_ptr>& tmp_pair = targets[j];
        TF&          tf2      = *(tmp_pair.first);
        coeffect_ptr cur_coef = tmp_pair.second;
        set(gene, tf1, tf2, cur_coef);
      }
    }   
  }
}

void ModifyingInteractions
::set(Gene& gene, TF& actor, TF& target, coeffect_ptr coef)
{
  site_ptr_vector& actorsites  = bindings->getSites(gene, actor);
  site_ptr_vector& targetsites = bindings->getSites(gene, target);

  gene_mods& gmods = *(mods[&gene]);
  vector<ModifyingInteraction>& v = gmods[&actor][&target];
  
  double max_dist     = coef->getMaxDistance();
  int    nactors      = actorsites.size();
  int    ntargets     = targetsites.size();

  for (int i=0; i<nactors; i++)
  {
    site_ptr actor_ptr  = actorsites[i];
    for (int j=0; j<ntargets; j++)
    {
      site_ptr target_ptr = targetsites[j];
      
      int m1 = actor_ptr->m;
      int n1 = actor_ptr->n;
      int m2 = target_ptr->m;  
      int n2 = target_ptr->n;  
      double d;
      
      int dm = abs(m1 - n2);
      int dn = abs(n1 - m2);
      
      if (dn <= dm)
        d = dn;
      else
        d = dm;
      
      bool  overlapped = (m1 < n2 && m2 < n1);
      double df        = coef->distFunc(d);
      if ( d < max_dist && !overlapped && df > 0)
      {
        ModifyingInteraction mod;
        mod.actor  = actor_ptr;
        mod.target = target_ptr;
        mod.distcoef = df;
        mod.coef = coef;
        
        v.push_back(mod);
      }
    }
  }
}

// [target][actor][nuc]
map<BindingSite*, map<TF*, vector<double> > > ModifyingInteractions::getEffq(Gene& gene)
{
  map<BindingSite*, map<TF*, vector<double> > > out;
  gene_mods& gmods = *(mods[&gene]);
  int ntfs   = tfs->size();
  for (int i = 0; i<ntfs; i++) // loop through actors
  {
    TF& tf1 = tfs->getTF(i);
    vector< pair<TF*, coeffect_ptr> >& targets = tf1.getTargets();
    int ntargets = targets.size();
    for (int j=0; j<ntargets; j++) // loop through targets
    {
      pair<TF*, coeffect_ptr>& tmp_pair = targets[j];
      TF&          tf2      = *(tmp_pair.first);
      coeffect_ptr cur_coef = tmp_pair.second;
      
      vector<ModifyingInteraction>& mod_vector = gmods[&tf1][&tf2];
      
      double efficiency = cur_coef->getEfficiency();
      int    coef_idx   = cur_coef->getIdx() - 1;
      
      int nmods = mod_vector.size();
      
      for (int k=0; k<nmods; k++)
      {
      
        ModifyingInteraction& mod = mod_vector[k];
        BindingSite& actor_site   = *mod.actor;
        BindingSite& target_site  = *mod.target;
        double distcoef           = mod.distcoef;
        
        vector<double>& actor_occupancy = actor_site.total_occupancy;
        vector<double>& start_occupancy = target_site.mode_occupancy[0];
        vector<double>& end_occupancy   = target_site.mode_occupancy[coef_idx];
        
        int n = actor_occupancy.size();

        if (out[&target_site][&tf1].size() == 0)
        {
          out[&target_site][&tf1].resize(n);
          for (int l=0; l<n; l++) 
            out[&target_site][&tf1][l] = 1.0;
        }
        
        for (int l=0; l<n; l++)
        {
          double reduction = actor_occupancy[l] * efficiency * distcoef;
          out[&target_site][&tf1][l] *= reduction;
          
        }
      }
    }
  }
  return out;
}
  
        
void ModifyingInteractions
::calc()
{
  initialize();
  int ngenes = genes->size();
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    int ntfs   = tfs->size();
    for (int i = 0; i<ntfs; i++) // loop through actors
    {
      TF& tf1 = tfs->getTF(i);
      vector< pair<TF*, coeffect_ptr> >& targets = tf1.getTargets();
      int ntargets = targets.size();
      for (int j=0; j<ntargets; j++)
      {
        pair<TF*, coeffect_ptr>& tmp_pair = targets[j];
        TF&          tf2      = *(tmp_pair.first);
        coeffect_ptr cur_coef = tmp_pair.second;
        calc(gene, tf1, tf2, cur_coef);
      }
    }   
  }
}

void ModifyingInteractions
::calc(Gene& gene)
{
  initialize(gene);

  int ntfs   = tfs->size();
  for (int i = 0; i<ntfs; i++) // loop through actors
  {
    TF& tf1 = tfs->getTF(i);
    vector< pair<TF*, coeffect_ptr> >& targets = tf1.getTargets();
    int ntargets = targets.size();
    for (int j=0; j<ntargets; j++)
    {
      pair<TF*, coeffect_ptr>& tmp_pair = targets[j];
      TF&          tf2      = *(tmp_pair.first);
      coeffect_ptr cur_coef = tmp_pair.second;
      calc(gene, tf1, tf2, cur_coef);
    }
  }   

}

void ModifyingInteractions
::calc(Gene& gene, TF& actor, TF& target, coeffect_ptr coef)
{
  gene_mods& gmods = *(mods[&gene]);
  vector<ModifyingInteraction>& mod_vector = gmods[&actor][&target];
  
  double efficiency = coef->getEfficiency();
  int    coef_idx   = coef->getIdx() - 1;
  
  int nmods = mod_vector.size();
  
  for (int i=0; i<nmods; i++)
  {

    ModifyingInteraction& mod = mod_vector[i];
    BindingSite& actor_site   = *mod.actor;
    BindingSite& target_site  = *mod.target;
    double distcoef           = mod.distcoef;
    
    vector<double>& actor_occupancy = actor_site.total_occupancy;
    vector<double>& start_occupancy = target_site.mode_occupancy[0];
    vector<double>& end_occupancy   = target_site.mode_occupancy[coef_idx];
    
    mod_f(actor_occupancy, start_occupancy, end_occupancy, efficiency, distcoef);
  } 
}

void ModifyingInteractions
::mod_f(vector<double>& actor_vec, 
        vector<double>& start_vec,
        vector<double>& end_vec,
        double ef, double d)
{
  int n = actor_vec.size();

  for (int i=0; i<n; i++)
  {
    double actor_occupancy = actor_vec[i];
    double reduction = actor_occupancy * ef * d;
    double current   = start_vec[i];
    double reduced   = current*(1-reduction);
    double increased = current - reduced;
    start_vec[i]     = reduced;
    end_vec[i]      += increased;
  }
}
    
void ModifyingInteractions
::initialize()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    initialize(gene);
  }
}

void ModifyingInteractions
::initialize(Gene& gene)
{
  int ntfs   = tfs->size();
  for (int j = 0; j<ntfs; j++) 
  {
    TF& tf = tfs->getTF(j);
    initialize(gene, tf);
  }
}
      

void ModifyingInteractions
::initialize(Gene& gene, TF& tf)
{
  int nmodes             = tf.getNumModes();
  site_ptr_vector& sites = bindings->getSites(gene, tf);

  int nsites = sites.size();

  for (int i=0; i<nsites; i++)
  {
    BindingSite* site = sites[i].get();
    vector<double>& total_occupancy = site->total_occupancy;
    site->mode_occupancy[0] = total_occupancy;
    int n = total_occupancy.size();
    for (int j=1; j<nmodes; j++)
    {
      vector<double>& tmp_occ = site->mode_occupancy[j];
      for (int k=0; k<n; k++)
        tmp_occ[k] = 0;
    }
  }
}




void ModifyingInteractions
::save(Gene& gene)
{

  gene_mods& gmods       = *(mods[&gene]);
  gene_mods& saved_gmods = *(saved_mods[&gene]);
  saved_gmods = gmods;

}

void ModifyingInteractions
::save()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    save(gene);
  }
}

void ModifyingInteractions
::clear(Gene& gene)
{
  gene_mods& gmods = *(mods[&gene]);
  gmods.clear();
}

void ModifyingInteractions
::clear()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    clear(gene);
  }
}

void ModifyingInteractions
::restore(Gene& gene)
{
  gene_mods& gmods       = *(mods[&gene]);
  gene_mods& saved_gmods = *(saved_mods[&gene]);
  gmods = saved_gmods;
}

void ModifyingInteractions
::restore()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    restore(gene);
  }
}


void ModifyingInteractions
::update()
{
  int ngenes = genes->size();
  for (int k=0; k<ngenes; k++)
  {
    Gene& gene = genes->getGene(k);
    update(gene);
  }
}

void ModifyingInteractions
::update(Gene& gene)
{
  gene_mods& gmods = *(mods[&gene]);
  int ntfs   = tfs->size();
  for (int i = 0; i<ntfs; i++) // loop through actors
  {
    TF& tf1 = tfs->getTF(i);
    vector< pair<TF*, coeffect_ptr> >& targets = tf1.getTargets();
    int ntargets = targets.size();
    for (int j=0; j<ntargets; j++)
    {
      pair<TF*, coeffect_ptr>& tmp_pair = targets[j];
      TF&          tf2      = *(tmp_pair.first);
      coeffect_ptr cur_coef = tmp_pair.second;
      gmods[&tf1][&tf2].clear();
      set(gene, tf1, tf2, cur_coef);
    }
  }   
}










