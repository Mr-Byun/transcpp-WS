/*********************************************************************************
*                                                                                *
*     subgroup.cpp                                                               *
*                                                                                *
*     Contains the subgroup classes and methods.                                 *
*                                                                                *
*********************************************************************************/

#include "subgroup.h"
#include <boost/foreach.hpp>
#include <algorithm>
#include <limits>
#include <cmath>
#include <climits>

#define foreach_ BOOST_FOREACH

// a fast implementation of log2 function, but with reduced precision...
inline float fast_log2 (float val)
{
   register int *const     exp_ptr = ((int*)&val);
   register int            x = *exp_ptr;
   register const int      log_2 = ((x >> 23) & 255) - 128;
   x &= ~(255 << 23);
   x += 127 << 23;
   *exp_ptr = x;

   val = ((-1.0f/3) * val + 2) * val - 2.0f/3;

   return ((val + log_2));
}

//inline float fast_log2 (float val)
//{
//   int * const    exp_ptr = reinterpret_cast <int * (&val);
//   int            x = *exp_ptr;
//   const int      log_2 = ((x  23) & 255) - 128;
//   x &= ~(255 << 23);
//   x += 127 << 23;
//   *exp_ptr = x;
//
//   val = ((-1.0f/3) * val + 2) * val - 2.0f/3;   // (1)
//
//   return (val + log_2);
//}
      
/********************************   Subgroup    *********************************/


/*
bool compareBindingSiteRight(BindingSite* a, BindingSite* b)
{
  return(a->n < b->n);
}

bool compareBindingSiteLeft(BindingSite* a, BindingSite* b)
{
  return(a->m > b->m);
}*/


/*    Constructors    */

Subgroup::Subgroup() {}

Subgroup::Subgroup(BindingSite* site, bindings_ptr b) 
{
  bindings  = b;
  addSite(site);
}


/*    Setters   */

void Subgroup::addSite(BindingSite* site)
{
  sites_f.push_back(site);
  if ( sites_f.size()==1 )
  {
    right_bound = site->n;
    left_bound  = site->m;
  } else
  {
    left_bound   = min(left_bound,  site->m);
    right_bound  = max(right_bound, site->n);
  }
}

void Subgroup::addSubgroup(Subgroup& sub)
{
  vector<BindingSite*>& subSites = sub.getSites();
  int nsites = subSites.size();
  for (int i=0; i<nsites; i++)
  {
    addSite( subSites[i] );
  }
}


/*    Getters   */

int Subgroup::size() const
{
  return sites_f.size();
}

vector<BindingSite*>& Subgroup::getSites()
{
  return sites_f;
}

int Subgroup::getRightBound() const { return right_bound; }

int Subgroup::getLeftBound() const { return left_bound; }


/*    Methods   */

bool Subgroup::overlaps(BindingSite& site)
{
  return ( site.m < right_bound && left_bound < site.n);
}

// check to see if this site coops with any sites already in group
bool Subgroup::checkCoop(BindingSite& site1)
{
  char o1,o2;
  TF& tf1 = *(site1.tf);
  
  if (!tf1.ncoops()) return false;
  
  int nsites = sites_f.size();
  int site1n = site1.n;
  int site1m = site1.m;  
  for (int i=0; i<nsites; i++)
  {
    BindingSite& site2 = *(sites_f[i]);
    int site2n = site2.n;
    int site2m = site2.m;  
    TF* tf2 = site2.tf;
    
    if (site1n <= site2m)
    {
      o1 = site1.orientation;
      o2 = site2.orientation;
    }
    else
    {
      o2 = site1.orientation;
      o1 = site2.orientation;
    }
    
    if (tf1.checkCoops(tf2, o1, o2))
    {
      coop_ptr cur_coop = tf1.getCoop(tf2);
      int dist1 = abs(site2m - site1n);
      int dist2 = abs(site1m - site2n);
      int dist  = min(dist1, dist2);
      double distfunc = cur_coop->distFunc(dist);
      if (distfunc > 0)
        return true;
    }
  }
  return false;
}

bool Subgroup::overlaps(BindingSite* site1, BindingSite* site2)
{
  return ( site1->m < site2->n && site2->m < site1->n);
}

void Subgroup::sort()
{
  int nsites = sites_f.size();
  
  std::sort(sites_f.begin(), sites_f.end(), compareBindingSiteRight);
  sites_r.resize(nsites);
  f2r.resize(nsites);
  for (int i=0; i<nsites; i++) 
    sites_r[i] = sites_f[nsites - i - 1];
  std::sort(sites_r.begin(), sites_r.end(), compareBindingSiteLeft);

  for (int i=0; i<nsites; i++)
  {
    for (int j=0; j<nsites; j++)
    {
      if (sites_f[i] == sites_r[j])
      {
        f2r[i] = j;
        continue;
      }
    }
  }
}

// for each site, find the last site that does not compete
void Subgroup::pre_process()
{
  int nsites = sites_f.size();
  int nnuc = bindings->getNnuc();
  
  ZF.resize(nsites+1);
  ZR.resize(nsites+1);
    
#ifdef VERYLARGENUMS
  ZF[0].log_Z.resize(nnuc);
  ZF[0].log_Zc.resize(nnuc);
  ZF[0].log_Znc.resize(nnuc);
  
  ZR[0].log_Z.resize(nnuc);
  ZR[0].log_Zc.resize(nnuc);
  ZR[0].log_Znc.resize(nnuc);
  for (int i=0; i<nnuc; i++)
  {
    ZF[0].log_Z[i] = 0;
    ZR[0].log_Z[i] = 0;
  }
#else  
  ZF[0].Z.resize(nnuc);
  ZF[0].Zc.resize(nnuc);
  ZF[0].Znc.resize(nnuc);
  
  ZR[0].Z.resize(nnuc);
  ZR[0].Zc.resize(nnuc);
  ZR[0].Znc.resize(nnuc);
  for (int i=0; i<nnuc; i++)
  {
    ZF[0].Z[i] = 1.0;
    ZR[0].Z[i] = 1.0;
  }
#endif
  
  for (int i=0; i<nsites; i++) // loop over all sites
  {
    int pindex = i + 1;
    
#ifdef VERYLARGENUMS
    ZF[pindex].log_Z.resize(nnuc);
    ZF[pindex].log_Zc.resize(nnuc, 0);
    ZF[pindex].log_Znc.resize(nnuc, 0);
    
    ZR[pindex].log_Z.resize(nnuc);
    ZR[pindex].log_Zc.resize(nnuc, 0);
    ZR[pindex].log_Znc.resize(nnuc, 0);
#else
    ZF[pindex].Z.resize(nnuc);
    ZF[pindex].Zc.resize(nnuc, 0);
    ZF[pindex].Znc.resize(nnuc, 0);
    
    ZR[pindex].Z.resize(nnuc);
    ZR[pindex].Zc.resize(nnuc, 0);
    ZR[pindex].Znc.resize(nnuc, 0);
#endif
    
    BindingSite* s1f = sites_f[i];
    BindingSite* s1r = sites_r[i];
    
    TF* tf1f = s1f->tf;
    TF* tf1r = s1r->tf;
    
    char o1f = s1f->orientation;
    char o1r = s1r->orientation;
    
    for (int j = 0; j<i; j++) // loop over previous sites
    {
      BindingSite* s2f = sites_f[j];
      BindingSite* s2r = sites_r[j];
      
      TF* tf2f = s2f->tf;
      TF* tf2r = s2r->tf;
      
      char o2f = s2f->orientation;
      char o2r = s2r->orientation;
      
      pre_process_pair(ZF, i, s1f, tf1f, o1f, j, s2f, tf2f, o2f); 
      pre_process_pair(ZR, i, s1r, tf1r, o1r, j, s2r, tf2r, o2r); 
    }
  }
}

void Subgroup::pre_process_pair(vector<Partition>& p, int i, BindingSite* s1, TF* tf1, char o1, int j, BindingSite* s2, TF* tf2, char o2)
{
  int pindex = i + 1;
  if (!overlaps(s1,s2))
  {
    p[pindex].last = j + 1;
    
    if (tf1->checkCoops(tf2, o2, o1))
    {
      coop_ptr cur_coop = tf1->getCoop(tf2);
      int dist1 = abs(s1->m - s2->n);
      int dist2 = abs(s2->m - s1->n);
      int dist = min(dist1,dist2);
      
      double distfunc = cur_coop->distFunc(dist);
      if (distfunc > 0)
      {
        p[pindex].coops = true;
        p[pindex].coop_site.push_back(j);
        p[pindex].coop_past.push_back(p[j+1].last);
        p[pindex].dist_coef.push_back(distfunc);
        p[pindex].coop.push_back(cur_coop);
      }
    }
  }
}
  

void Subgroup::
iterate_partition(vector<Partition>& p, vector<BindingSite*>& sites, int site_index)
{
  int pindex           = site_index + 1;

  Partition& cur_part  = p[pindex];
  Partition& last_part = p[cur_part.last];
  Partition& init_part = p[site_index];
    
#ifdef VERYLARGENUMS
  vector<double>& cur_part_logZ      = cur_part.log_Z;
  vector<double>& cur_part_logZnc    = cur_part.log_Znc;
  vector<double>& cur_part_logZc     = cur_part.log_Zc;
#elif defined LARGENUMS
  vector<long double>& cur_part_Z   = cur_part.Z;
  vector<long double>& cur_part_Znc = cur_part.Znc;
  vector<long double>& cur_part_Zc  = cur_part.Zc;
#else
  vector<double>& cur_part_Z   = cur_part.Z;
  vector<double>& cur_part_Znc = cur_part.Znc;
  vector<double>& cur_part_Zc  = cur_part.Zc;
#endif

  int ncoops = cur_part.coop_site.size();
    
  vector<double>& kv   = sites[site_index]->kv;
  
  int nnuc = bindings->getNnuc();
  
#ifdef VERYLARGENUMS
  for (int i=0; i<nnuc; i++)
  {
    double cur_kv = kv[i];
    double init_logZ = init_part.log_Z[i];
    double last_logZ = last_part.log_Z[i];
    double new_Z_ratio   = 1;
    double new_Znc_ratio = exp2(last_logZ - init_logZ)*cur_kv;
    double new_Zc_ratio  = 0;
    
    for (int j=0; j<ncoops; j++)
    {
      double kcoop      = cur_part.coop[j]->getK();
      int coop_site = cur_part.coop_site[j]; // the site it coops with
      int coop_past = cur_part.coop_past[j];
      double dfunk      = cur_part.dist_coef[j];
      kcoop *= dfunk;
      
      double cur_coopkv = sites[coop_site]->kv[i]; // kv of this coop site
      last_logZ         = p[coop_past].log_Z[i];   // last z before coop
      new_Zc_ratio += exp2(last_logZ - init_logZ)*cur_coopkv*cur_kv*kcoop;
    }
    
    new_Z_ratio = 1 + new_Znc_ratio + new_Zc_ratio;
    cur_part_logZ[i]   = log2(new_Z_ratio)  + init_logZ;
    cur_part_logZnc[i] = log2(new_Znc_ratio) + init_logZ;
    cur_part_logZc[i]  = log2(new_Zc_ratio)  + init_logZ;
  }
#elif defined LARGENUMS
  for (int i=0; i<nnuc; i++)
  {
    long double cur_kv = kv[i];
    long double init_Z = init_part.Z[i];
    long double last_Z = last_part.Z[i];
    long double new_Z  = last_Z*cur_kv;
    
    cur_part_Z[i]   = init_Z + new_Z;
    cur_part_Znc[i] = new_Z;
    cur_part_Zc[i]  = 0;
  }
    
  for (int j=0; j<ncoops; j++)
  {
    long double kcoop      = cur_part.coop[j]->getK();
    unsigned int coop_site = cur_part.coop_site[j]; // the site it coops with
    unsigned int coop_past = cur_part.coop_past[j];
    long double dfunk      = cur_part.dist_coef[j];
    kcoop *= dfunk;
    
    vector<double>& coopkv = sites[coop_site]->kv;
    
    for (int k=0; k<nnuc; k++)
    {
      long double cur_coopkv = coopkv[k];
      long double cur_kv     = kv[k];
      long double last_Z     = p[coop_past].Z[k];
      long double weight     = last_Z*cur_coopkv*cur_kv*kcoop;
      cur_part_Z[k]  += weight;
      cur_part_Zc[k] += weight;
    }
  }
#else
  for (int i=0; i<nnuc; i++)
  {
    double cur_kv = kv[i];
    double init_Z = init_part.Z[i];
    double last_Z = last_part.Z[i];
    double new_Z  = last_Z*cur_kv;
    
    cur_part_Z[i]   = init_Z + new_Z;
    cur_part_Znc[i] = new_Z;
    cur_part_Zc[i]  = 0;
  }
    
  for (int j=0; j<ncoops; j++)
  {
    double kcoop      = cur_part.coop[j]->getK();
    int coop_site = cur_part.coop_site[j]; // the site it coops with
    int coop_past = cur_part.coop_past[j];
    double dfunk      = cur_part.dist_coef[j];
    kcoop *= dfunk;
    
    vector<double>& coopkv = sites[coop_site]->kv;
    
    for (int k=0; k<nnuc; k++)
    {
      double cur_coopkv = coopkv[k];
      double cur_kv     = kv[k];
      double last_Z     = p[coop_past].Z[k];
      double weight     = last_Z*cur_coopkv*cur_kv*kcoop;
      cur_part_Z[k]  += weight;
      cur_part_Zc[k] += weight;
    }
  }
#endif
}

     
void Subgroup::occupancy()
{
  int nsites = sites_f.size();
  int nnuc   = bindings->getNnuc();
  
  for (int i=0; i<nsites; i++)
  {
    iterate_partition(ZF, sites_f, i); // Z0 now holds all the partitions
    iterate_partition(ZR, sites_r, i);
  }
   
#ifdef VERYLARGENUMS
  vector<double>& log_Z = ZF[nsites].log_Z;

  for (int i=0; i<nsites; i++)
  { 
    int r_idx = f2r[i];
    
    vector<double>& log_zfnc = ZF[i+1].log_Znc;
    vector<double>& log_zrnc = ZR[r_idx+1].log_Znc;
    
    vector<double>& log_zfc = ZF[i+1].log_Zc;
    vector<double>& log_zrc = ZR[r_idx+1].log_Zc;
      
    BindingSite* site = sites_f[i];
    
    for (int j=0; j<nnuc; j++)
    {
      
      double kv = site->kv[j];
      if ( kv == 0)
      {
        site->total_occupancy[j]   = 0;
        site->mode_occupancy[0][j] = 0;
      }
      else
      {  
        double nc = exp2(log_zfnc[j] + log_zrnc[j] - log_Z[j])/kv;
        double cr = exp2(log_zfnc[j] + log_zrc[j]  - log_Z[j])/kv;
        double cf = exp2(log_zfc[j]  + log_zrnc[j] - log_Z[j])/kv;
        
        double f  = nc + cr + cf;
        
        if (std::isnan(f)) error("f is NaN");
        site->total_occupancy[j]   = f;
        site->mode_occupancy[0][j] = f;
      }
    }
  }
#elif defined LARGENUMS
  vector<long double>& Z = ZF[nsites].Z;

  for (int i=0; i<nsites; i++)
  { 
    int r_idx = f2r[i];
    
    vector<long double>& zfnc = ZF[i+1].Znc;
    vector<long double>& zrnc = ZR[r_idx+1].Znc;
    
    vector<long double>& zfc = ZF[i+1].Zc;
    vector<long double>& zrc = ZR[r_idx+1].Zc;
      
    BindingSite* site = sites_f[i];
    
    for (int j=0; j<nnuc; j++)
    {
      
      long double kv = site->kv[j];
      if ( kv == 0)
      {
        site->total_occupancy[j]   = 0;
        site->mode_occupancy[0][j] = 0;
      }
      else
      {  
        long double denom = (zfnc[j] + zfc[j])*(zrnc[j] + zrc[j]) - zrc[j]*zfc[j];
        double f  = denom/(Z[j]*kv);
        
        if (std::isnan(f)) error("The partition function overflowed 'long double'. Recompile with VERYLARGENUMS=ON");
        site->total_occupancy[j]   = f;
        site->mode_occupancy[0][j] = f;
      }
    }
  }
#else
  vector<double>& Z = ZF[nsites].Z;

  for (int i=0; i<nsites; i++)
  { 
    int r_idx = f2r[i];
    
    vector<double>& zfnc = ZF[i+1].Znc;
    vector<double>& zrnc = ZR[r_idx+1].Znc;
    
    vector<double>& zfc = ZF[i+1].Zc;
    vector<double>& zrc = ZR[r_idx+1].Zc;
      
    BindingSite* site = sites_f[i];
    
    for (int j=0; j<nnuc; j++)
    {
      
      double kv = site->kv[j];
      if ( kv == 0)
      {
        site->total_occupancy[j]   = 0;
        site->mode_occupancy[0][j] = 0;
      }
      else
      {  
        double denom = (zfnc[j] + zfc[j])*(zrnc[j] + zrc[j]) - zrc[j]*zfc[j];
        double f  = denom/(Z[j]*kv);
        
        if (std::isnan(f)) error("The partition function overflowed 'long double'. Recompile with LARGENUMS=ON");
        site->total_occupancy[j]   = f;
        site->mode_occupancy[0][j] = f;
      }
    }
  }
#endif
}
      

/*    Print   */

void Subgroup::print(ostream& os)
{
  os << setw(10) << "f_index";
  os << setw(10) << "r_index";
  
  printSiteHeader(os);
  
  os << setw(10) << "last_f" << setw(10) << "last_r" << endl;
  
  int nsites = sites_f.size();
  cerr << nsites << endl;
  for (int i=0; i<nsites; i++)
  {
    os << setw(10) << i;
    os << setw(10) << f2r[i];
    
    printSite(*sites_f[i],os);
    os << setprecision(3) << setw(10);
    os << ZF[i+1].last;
    os << setw(10) << ZR[f2r[i]+1].last << endl;
    os << endl;
  }
}


/******************************   Subgroups   **********************************/

Subgroups::Subgroups() {}

void Subgroups::clear()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
    clear(genes->getGene(i));
}

void Subgroups::clear(Gene& gene)
{
  list<Subgroup>& ggroups = *(groups[&gene]);
  ggroups.clear();
}


void Subgroups::create(genes_ptr g, tfs_ptr t, bindings_ptr b, mode_ptr m)
{
  genes     = g;
  tfs       = t;
  bindings  = b;
  mode      = m;
  
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    gene_groups_ptr ggroups(new list<Subgroup>);
    gene_groups_ptr saved_ggroups(new list<Subgroup>);
    
    groups[&gene]       = ggroups;
    saved_groups[&gene] = saved_ggroups;
  }
  
  update();
}

void Subgroups::addSites(Gene& gene)
{
  map<TF*, site_ptr_vector>& gene_sites = bindings->getSites(gene);
  list<Subgroup>& gene_groups = *(groups[&gene]);
  
  int ntfs = tfs->size();
  for(int i=0; i<ntfs; i++)
  {  
    TF& tf = tfs->getTF(i);
    site_ptr_vector& sites = gene_sites[&tf];
    int nsites = sites.size();
    for (int j=0; j<nsites; j++)
      addSite(gene_groups, sites[j]);
  }
  
  list<Subgroup>::iterator i;
  for (i=gene_groups.begin(); i != gene_groups.end(); ++i)
  {
    i->sort();
    i->pre_process();
  }
}

void Subgroups::addSite(list<Subgroup>& gene_groups, site_ptr site)
{
  list<Subgroup>::iterator i;  // iterate through list
  list<Subgroup>::iterator j;  // point to group that site was added to
   
  bool added = false;  
  i = gene_groups.begin();
  while (i != gene_groups.end())
  {
    BindingSite& site_ref = *site;
    
    /* I check every subgroup. I add the site to the first one that it overlaps or
    coops with. If it coops with additional ones I merge them. */
    bool condition = i->overlaps(site_ref) || i->checkCoop(site_ref);
    
    if (condition && added==false)
    {
      i->addSite(site.get());
      added = true;
      j = i;
      i++;
    } else if ( condition && added==true)
    {
      j->addSubgroup(*i);
      i = gene_groups.erase(i);
    } else {
      i++;
    }
  }
  
  if (added==false)
    gene_groups.push_back(Subgroup(site.get(), bindings));
}
    
  

void Subgroups::update()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    update(gene);
  }
}

void Subgroups::update(Gene& gene)
{
  list<Subgroup>& ggroups = *(groups[&gene]);
  ggroups.clear();
  addSites(gene);
}
 

void Subgroups::calc_f()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
    calc_f(genes->getGene(i));
}

void Subgroups::calc_f(Gene& gene)
{
  list<Subgroup>& gene_groups = *(groups[&gene]);
  list<Subgroup>::iterator i;
  for (i=gene_groups.begin(); i != gene_groups.end(); ++i)
    i->occupancy();
}


void Subgroups::save()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
    save(genes->getGene(i));
}

void Subgroups::save(Gene& gene)
{
  list<Subgroup>& gene_groups       = *(groups[&gene]);
  list<Subgroup>& saved_gene_groups = *(saved_groups[&gene]);
  saved_gene_groups = gene_groups;
}

void Subgroups::restore()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
    restore(genes->getGene(i));
}

void Subgroups::restore(Gene& gene)
{
  list<Subgroup>& gene_groups       = *(groups[&gene]);
  list<Subgroup>& saved_gene_groups = *(saved_groups[&gene]);
  gene_groups = saved_gene_groups;
}

void Subgroups::print(ostream& os)
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    os << gene.getName();
    print(gene, os);
  }
}

void Subgroups::print(Gene& gene, ostream& os)
{
  list<Subgroup>& gene_groups = *(groups[&gene]);
  int counter = 1;
  list<Subgroup>::iterator i;
  
  for (i=gene_groups.begin(); i != gene_groups.end(); ++i)
  {
    os << "Subroup " << counter++ << endl;
    os << "Left bound:  " << i->getLeftBound()  << endl;
    os << "Right bound: " << i->getRightBound() << endl;
    i->print(os);
  }
}
