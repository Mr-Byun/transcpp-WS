/*********************************************************************************
*                                                                                *
*     nuclei.cpp                                                                 *
*                                                                                *
*     A nuclei object contains all nuclei which share the same TFs. This means   *
*     they share subgroups and quenching interactions.                           *
*                                                                                *
*     In the original code, many of the calculations on occupancy and quenching  *
*     were multiplications by 0. This class was created to prevent such needless *
*     calculations.                                                              *
*                                                                                *
*********************************************************************************/

#include "nuclei.h"

#include <boost/foreach.hpp>
#include <limits>

# define foreach_ BOOST_FOREACH

/*    Constructors    */

Nuclei::Nuclei() :
  bindings(bindings_ptr(new Bindings)),
  subgroups(subgroups_ptr(new Subgroups)),
  quenching(quenching_ptr(new QuenchingInteractions)),
  coeffects(modifying_ptr(new ModifyingInteractions))
{}


Nuclei::Nuclei(Organism* parent, tfs_ptr t) :
  bindings(bindings_ptr(new Bindings)),
  subgroups(subgroups_ptr(new Subgroups)),
  quenching(quenching_ptr(new QuenchingInteractions)),
  coeffects(modifying_ptr(new ModifyingInteractions))
{
  n           = 0;
  tfs         = t;
  genes       = parent->getGenes();
  tfdata      = parent->getTFData();
  distances   = parent->getDistances();
  promoters   = parent->getPromoters();
  mode        = parent->getMode();
  competition = parent->getCompetition();
  chromatin   = parent->getChromatin();
  
  competition_mode = mode->getCompetition();
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    penalty[&gene] = 1.0;
  }
}

void Nuclei::setParent(Organism* parent)
{
  n           = 0;
  tfs         = parent->getTFs();
  genes       = parent->getGenes();
  tfdata      = parent->getTFData();
  distances   = parent->getDistances();
  promoters   = parent->getPromoters();
  mode        = parent->getMode();
  competition = parent->getCompetition();
  chromatin   = parent->getChromatin();
  
  competition_mode = mode->getCompetition();
}
  

/*    Setters    */

void Nuclei::clear()
{
  n = 0;
  IDs.clear();
  Ns.clear();
  Rs.clear();
  competition_map.clear();
}

void Nuclei::setTFs(tfs_ptr x)             {tfs       = x;}
void Nuclei::setGenes(genes_ptr x)         {genes     = x;}
void Nuclei::setBindings(bindings_ptr x)   {bindings  = x;}
void Nuclei::setTFData(table_ptr x)        {tfdata    = x;}
void Nuclei::setDistances(distances_ptr x) {distances = x;}
void Nuclei::setPromoters(promoters_ptr x) {promoters = x;}

void Nuclei::create(ptree& pt)
{
  createBindings(pt);
  createSubgroups();
  createCoeffects();
  createQuenching();

  calcOccupancy();
  calcCoeffects();
  calcQuenching();
  calcR();
}

void Nuclei::create()
{
  createBindings();
  createSubgroups();
  createCoeffects();
  createQuenching();

  calcOccupancy();
  calcCoeffects();
  calcQuenching();
  calcR();
}


void Nuclei::createBindings(ptree& pt)
{
  bindings->setGenes(genes);
  bindings->setTFs(tfs);
  bindings->setTFData(tfdata); 
  bindings->setMode(mode); 
  bindings->setChromatin(chromatin);
  
  int nids = IDs.size();
  for (int i=0; i<nids; i++)
    bindings->addNuc(IDs[i]);

  if(mode->getBindingSiteList())
    bindings->createSites(pt);
  else
    {
      bindings->createScores();
      bindings->createSites();
    }
}

void Nuclei::createBindings()
{
  bindings->setGenes(genes);
  bindings->setTFs(tfs);
  bindings->setTFData(tfdata); 
  bindings->setMode(mode); 
  bindings->setChromatin(chromatin);
  
  int nids = IDs.size();
  for (int i=0; i<nids; i++)
    bindings->addNuc(IDs[i]);

  bindings->createScores();
  bindings->createSites();
}

void Nuclei::createSubgroups()
{
  subgroups->create(genes, tfs, bindings, mode);
}

void Nuclei::createCoeffects()
{
  coeffects->create(genes, tfs, bindings, distances);
}

void Nuclei::createQuenching()
{
  quenching->create(genes, tfs, bindings, distances);
}
  
void Nuclei::addNuc(string id)
{
  IDs.push_back(id);
  n=IDs.size();
  
  int ngenes = genes->size();
  
  if (!competition_mode)
  {
    for (int i=0; i<ngenes; i++)
    {
      Gene& gene = genes->getGene(i);
      Ns[&gene].resize(n);
      Rs[&gene].resize(n);
    }
  }
  else
  {
    int window = competition->getWindow();
    int shift  = competition->getShift();
    
    for (int i=0; i<ngenes; i++)
    {
      Gene& gene = genes->getGene(i);
      Rs[&gene].resize(n);
      int length = gene.length() + window - shift;
      competition_data& gcomp = competition_map[&gene];
      gcomp.nwindows = max(1, length/shift + (length % shift != 0));
      //cerr << "new window size is " << window << endl;
      //gcomp.nwindows = length/shift + 2*window/shift - 2;
      
      gcomp.N_2D.resize(gcomp.nwindows);
      gcomp.T_2D.resize(gcomp.nwindows);
      gcomp.R_2D.resize(gcomp.nwindows);
      
      for (int j=0; j<gcomp.nwindows; j++)
      {
        gcomp.N_2D[j].resize(n);
        gcomp.T_2D[j].resize(n);
        gcomp.R_2D[j].resize(n);
      }
      
      gcomp.delta_N.resize(n);
      gcomp.total_N.resize(n);
    }
  }   
}

void Nuclei::resizeWindows()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    resizeWindow(gene);
  }
}

void Nuclei::resizeWindow(Gene& gene)
{
  int window = competition->getWindow();
  int shift  = competition->getShift();

  int length = gene.length() + window - shift;
  competition_data& gcomp = competition_map[&gene];
  gcomp.nwindows = max(1, length/shift + (length % shift != 0));
  //gcomp.nwindows = length/shift + 2*window/shift - 2;
  
  gcomp.N_2D.resize(gcomp.nwindows);
  gcomp.T_2D.resize(gcomp.nwindows);
  gcomp.R_2D.resize(gcomp.nwindows);
  
  for (int j=0; j<gcomp.nwindows; j++)
  {
    gcomp.N_2D[j].resize(n);
    gcomp.T_2D[j].resize(n);
    gcomp.R_2D[j].resize(n);
  }
}

/*    Getters   */

tfs_ptr       Nuclei::getTFs()          {return tfs;      }
genes_ptr     Nuclei::getGenes()        {return genes;    }
bindings_ptr  Nuclei::getBindings()     {return bindings; }
table_ptr     Nuclei::getTFData()       {return tfdata;   }
distances_ptr Nuclei::getDistances()    {return distances;}
promoters_ptr Nuclei::getPromoters()    {return promoters;}
quenching_ptr Nuclei::getQuenching()    {return quenching;}
subgroups_ptr Nuclei::getSubgroups()    {return subgroups;}
  
double& Nuclei::getRate(Gene& gene, string& id)
{
  int nids = IDs.size();
  for (int i=0; i<nids; i++)
  {
    if (id == IDs[i])
    {
      return Rs[&gene][i];
    }
  }
  stringstream err;
  err << "ERROR: could not find id in this set of nuclei!" << endl;
  error(err.str());
  return Rs[&gene][0]; // you will never get here!
}

/*    Methods   */

bool Nuclei::compareTFs(tfs_ptr t)
{
  if (tfs->size() != t->size())
    return false;
  else
  {
    for (int i=0; i<tfs->size(); i++)
    {
      if (tfs->getTFptr(i) != t->getTFptr(i))
        return false;
    }
  }
  return true;
}

/*
void Nuclei::calcN()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    calcN(gene);
  }
}*/


void Nuclei::calcN(Gene& gene)
{
  vector<double>& tNs = Ns[&gene];
  tNs.resize(n);
  for (int i=0; i<n;i++)
    tNs[i]=0;
  
  int ntfs = tfs->size();
  
  for (int i = 0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    
    if (tf.neverActivates()) continue;
    
    vector<double> coefs    = tf.getCoefs();
    site_ptr_vector& tsites = bindings->getSites(gene,tf);
    int ntsites = tsites.size();
    int nmodes  = coefs.size();
    for (int j=0; j<nmodes; j++)
    {
      double efficiency = coefs[j];
      if (efficiency <= 0) continue;

      for (int k=0; k<ntsites; k++)
      {
        vector<double>& eff_occ = tsites[k]->effective_occupancy[j];
        for (int l=0; l<n; l++)
        {
          tNs[l] += eff_occ[l]*efficiency;
        }
      }
    }
  }
}

void Nuclei::calcR()
{
  int ngenes = genes->size();
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes->getGene(i);
    calcR(gene);
  }
}

void Nuclei::calcR(Gene& gene)
{
  vector<double>& tNs = Ns[&gene];
  vector<double>& tRs = Rs[&gene];
  penalty[&gene] = 1.0;
  
  if (!competition_mode)
  {
    calcN(gene);
    for (int i=0; i<n; i++)
      tRs[i] = gene.getRate(tNs[i]);
  }
  else
    calcR2(gene);
}

        
/* THIS FUNCTION WILL ONLY WORK IF EFS TAKE ON VALUES BETWEEN 
0 AND 1. IF MORE ACTIVATION IS NEEDED USE Q INSTEAD */
void Nuclei::getRDist(Gene& gene, vector<BindingSite*>& sites, int start_idx, int end_idx, vector<double>& out)
{
  // make a convenient matrix for holding occupanies and efficiencies
  // p indexes [site][nuc]
  vector<vector<double> > p;
  vector<double>          effs;
  //cerr << "1" << endl;
  for (int i = start_idx; i <= end_idx; i++)
  {
    BindingSite& site    = *sites[i];
    TF& tf               = *site.tf;

    vector<double> coefs = tf.getCoefs();        
    int nmodes           = coefs.size();
    for (int j=0; j<nmodes; j++)
    {
      double efficiency = coefs[j];
      if (efficiency <= 0) continue;
      p.push_back(site.effective_occupancy[j]);
      int t = p.size() - 1;
      for (int k=0; k<n; k++)
      {
        p[t][k] *= efficiency;
        if (p[t][k] < 0) cerr << "efficiency less than 0!" << endl;
        if (p[t][k] > 1) cerr << "efficiency greater than 1!" << endl;
      }
    }
  }
  
  //cerr << "2" << endl;
  
  vector<vector<double> > n_dist;
  n_dist.resize(p.size()+1);
  for (int i=0; i< (int) n_dist.size(); i++)
    n_dist[i].resize(n);
  
  calc_pascal_2D(p, n_dist); 
  //cerr << "3" << endl;
  for (int j=0; j<n; j++)
  {
    out[j] = 0;
    double sum = 0;
    for (int i=0; i< (int) n_dist.size(); i++)
    {
      //cerr << n_dist[i][j] << endl;
      out[j] += n_dist[i][j]*gene.getRate(i);
      sum += n_dist[i][j];
    }
    //cerr << "i:   " << j << endl;
    //cerr << "sum: " << sum << endl;
    //cerr << "out: " << out[j] << endl;
  }

}

void Nuclei::calc_pascal(vector<double> in)
{
  int length = in.size() + 1;
  vector<vector<double> > out;
  out.resize(length);
  for (int i=0; i<length; i++)
    out[i].resize(i+1);

  out[0][0] = 1;
  for (int i=1; i<length; i++)
  {
    vector<double>& outi  = out[i];
    vector<double>& outi1 = out[i-1];
    double p = in[i-1];
    for (int j=0; j<i; j++)
      outi[j] += outi1[j]*(1-p);
    for (int j=1; j<=i; j++)
      outi[j] += outi1[j-1]*p;
  }
  
  for (int i=0; i<length; i++)
  {
    
    for (int j=0; j<=i; j++)
    {
      cerr << out[i][j] << " ";
    }
    cerr << endl;
  }
}

// in is [site#][nuc]
void Nuclei::calc_pascal_2D(vector<vector<double> >& in, vector<vector<double> >& ret)
{
  int length = in.size() + 1;
  // out is [site#][nuc][row]
  vector<vector<vector<double> > > out;
  out.resize(length);
  for (int i=0; i<length; i++)
  {
    out[i].resize(i+1);
    for (int j=0; j<=i; j++)
      out[i][j].resize(n);
  }

  for (int i=0; i<n; i++)
    out[0][0][i] = 1;

  for (int i=1; i<length; i++)
  {
    for (int k=0; k<n; k++)
    {
      vector<vector<double> >& outi  = out[i];
      vector<vector<double> >& outi1 = out[i-1];
      double p = in[i-1][k];
      for (int j=0; j<i; j++)
        outi[j][k] += outi1[j][k]*(1-p);
      for (int j=1; j<=i; j++)
        outi[j][k] += outi1[j-1][k]*p;
    }
  }
  
  for (int k=0; k<n; k++)
  {
    for (int i=0; i<length; i++)
    {
      ret[i][k] = out[length-1][i][k];
      //cerr << ret[i][k] << " ";
    }
    //cerr << endl;
  }
}
    
  

void Nuclei::calcR2(Gene& gene)
{
  competition_data& gcomp = competition_map[&gene];
  
  vector<BindingSite*>& sites_f = bindings->getFsites(gene);
  vector<BindingSite*>& sites_r = bindings->getRsites(gene);
 
  int window         = competition->getWindow();
  int shift          = competition->getShift();
  double specificity = competition->getSpecificity();
  double threshold   = competition->getThreshold();
  double background  = competition->getBackground();
  bool product       = competition->getProduct();
  double S           = competition->getS();
  double interaction_strength = competition->getInteractionStrength();
  
  int nsites = sites_f.size();
  
  int prime5 = 0 - window;
  int prime3 = prime5 + window;
  
  int site_f_idx = 0;
  int site_r_idx = nsites - 1;
  
  for (int j=0; j<n; j++)
  {
    gcomp.delta_N[j] = 0.0;
    gcomp.total_N[j] = background;
  }
    
  int nwindows = gcomp.nwindows;
  
  for (int i=0; i<nwindows; i++)
  {
    vector<double>& sub_N = gcomp.N_2D[i];
    vector<double>& sub_R = gcomp.R_2D[i];
    vector<double>& sub_T = gcomp.T_2D[i];
    
    prime3 += shift;
    prime5 += shift;
    
    // see what new sites have been added in this window
    bool in_window = true;
    
    while (in_window && site_f_idx < nsites)
    {
      BindingSite& site = *sites_f[site_f_idx];
      int pos = (site.m+site.n)/2; 
      if (pos < prime3)
      {
        TF& tf = *site.tf;
        vector<double> coefs    = tf.getCoefs();                
        int nmodes  = coefs.size();
        
        for (int j=0; j<nmodes; j++)
        {
          double efficiency = coefs[j];
          
          if (efficiency <= 0) continue;
          
          vector<double>& eff_occ = site.effective_occupancy[j];
          
          for (int k=0; k<n; k++)
            gcomp.delta_N[k] +=  eff_occ[k]*efficiency;
        }
        site_f_idx++;
      }
      else
        in_window = false;
    }
    
    // see what sites have been removed
    bool out_window = true;
    
    while (out_window && site_r_idx >= 0)
    {
      BindingSite& site = *sites_r[site_r_idx];
      int pos = (site.m+site.n)/2; 
      if (pos < prime5)
      {
        TF& tf = *site.tf;
        vector<double> coefs    = tf.getCoefs();                
        int nmodes  = coefs.size();
        
        for (int j=0; j<nmodes; j++)
        {
          double efficiency = coefs[j];
          
          if (efficiency <= 0) continue;
          
          vector<double>& eff_occ = site.effective_occupancy[j];
          
          for (int k=0; k<n; k++)
            gcomp.delta_N[k] -=  eff_occ[k]*efficiency;
        }
        site_r_idx--;
      }
      else
        out_window = false;
    }
   
    //cerr << "site_f_idx: " << site_f_idx << endl;
    //cerr << "site_r_idx: " << site_r_idx << endl;
    
    //BindingSite& sitef = *sites_f[site_f_idx-1];
    //BindingSite& sitel = *sites_r[site_r_idx];
    //int end_idx = sitef.index_in_ordered_f;
    //int start_idx   = sitel.index_in_ordered_f;

    //cerr << "gene:  " << gene.getName() << endl;
    //cerr << "start: " << start_idx << " of " << sites_f.size() << endl;
    //cerr << "end:   " << end_idx   << " of " << sites_f.size() << endl;
    //getRDist(gene, sites_f, start_idx, end_idx, sub_R);
    
    if (product)
    {
      
      for (int j=0; j<n;j++)
      {
        double& nj = sub_N[j];
        nj = gcomp.delta_N[j];
        sub_R[j] = gene.getRate(nj);
        double p = pow(S, nj);
        sub_T[j] = p*interaction_strength;
        gcomp.total_N[j] += p*interaction_strength;
      }
    }
    else
    {
      for (int j=0; j<n;j++)
      {
        double& nj = sub_N[j];
        nj = gcomp.delta_N[j];
        sub_R[j] = gene.getRate(nj);
  
        if ( nj >= threshold )
        {
          double p = pow(nj, specificity);
          sub_T[j] = p*interaction_strength;
          gcomp.total_N[j] += p*interaction_strength;
        }
      }
    }
  }
  
  vector<double>& tRs = Rs[&gene];
  for (int j=0; j<n;j++)
    tRs[j] = 0.0;
  
  for (int i=0; i<nwindows; i++)
  {
    vector<double>& sub_N = gcomp.N_2D[i];
    vector<double>& sub_R = gcomp.R_2D[i];
    vector<double>& sub_T = gcomp.T_2D[i];
    
    if (product)
    {
      for (int j=0; j<n;j++)
      {
        double nj = sub_N[j];
        if ( nj < threshold )
        {
          sub_T[j] = 0.0;
          tRs[j] += 0.0;
        }
        else
        {
          sub_T[j] /= gcomp.total_N[j];
          tRs[j] += sub_R[j] * sub_T[j];
        }
      }
    }
    else
    {
      for (int j=0; j<n;j++)
      {
        double nj = sub_N[j];
        if ( nj < threshold )
        {
          sub_T[j] = 0.0;
          tRs[j] += 0.0;
        }
        else
        {
          sub_T[j] /= gcomp.total_N[j];
          tRs[j] += sub_R[j] * sub_T[j];
        }
      }
    }
  }
  
  //double penalty_threshold     = 0.5;
  //
  //double enhancer_lower_size   = 500;
  //double enhancer_higher_size  = 2000;
  //double enhancer_size_penalty = 1;
  //enhancer_higher_size /= shift;
  //enhancer_lower_size  /= shift;
  //
  //double enhancer_num_lower    = 3;
  //double enhancer_num_higher   = 6;
  //double enhancer_num_penalty  = 1;
  //
  //vector<int> penalty_windows; // number of windows driving more than 0.5
  //penalty_windows.resize(n);
  //vector<int> penalty_nuc;     // number of nuclei driving more than 0.5
  //vector<bool> drives;
  //drives.resize(nwindows);
  //penalty_nuc.resize(nwindows);
  
  // get the max rate
  //double max_rate = 0;
  //for (int i=0; i<n; i++)
  //{
  //  double rate = tRs[i];
  //  if (rate > max_rate)
  //    max_rate = rate;
  //}
  //
  //double penalty_cutoff = max_rate * penalty_threshold;
  //for (int i = 0; i<nwindows; i++)
  //{
  //  drives[i] = 0;
  //  vector<double>& sub_R = gcomp.R_2D[i];
  //  for (int j=0; j<n;j++)
  //  {
  //    double rj = sub_R[j];
  //    if (rj > penalty_cutoff)
  //    {
  //      drives[i] = 1;
  //      //penalty_windows[j]++;
  //      //penalty_nuc[i]++;
  //    } 
  //  }
  //}
  
  // where driven
  /*for (int i = 0; i<nwindows; i++)
  {
    if (!drives[i]) continue;
    double pn = penalty_nuc[i];
    //cerr << "penalty_nuc[" << i << "] = " << pn << endl;
    if (pn > enhancer_num_higher)
    {
      enhancer_num_penalty *= pn / enhancer_num_higher;
      cerr << "at window " << i << " number of nuc " << pn << " is higher than " << enhancer_num_higher << " for penalty " << pn / enhancer_num_higher << endl;
    }
    else if (pn < enhancer_num_lower)
    {
      enhancer_num_penalty *= enhancer_num_lower / pn;
      cerr << "at window " << i << " number of nuc " << pn << " is lower than " << enhancer_num_lower << " for penalty " << enhancer_num_lower / pn << endl;
    }
  }*/
  
  // where expresses, should not in more than n windows
  /*for (int i = 0; i<n; i++)
  {
    if ( tRs[i] <= penalty_cutoff) continue;
    double pw = penalty_windows[i];
    //cerr << "penalty_windows[" << i << "] = " << pw << endl;
    if (pw > enhancer_higher_size)
    {
      enhancer_size_penalty *= pw / enhancer_higher_size;
      cerr << "at nuc " << n << " number of windows " << pw << " is lower than " << enhancer_higher_size << " for penalty " << pw / enhancer_higher_size << endl;
    }
    else if (pw < enhancer_lower_size)
    {
      enhancer_size_penalty *= enhancer_lower_size / pw;
      cerr << "at nuc " << n << " number of windows " << pw <<  "is lower than " << enhancer_lower_size << " for penalty " << enhancer_lower_size / pw << endl;
    }
    
    cerr << enhancer_size_penalty << endl;
  }*/
  
  //cerr << "enhancer_size_penalty: " << enhancer_size_penalty << endl;
  //cerr << "enhancer_num_penalty:  " << enhancer_num_penalty << endl;
  
  //int total_drives = 0;
  //for (int i=0; i<nwindows; i++)
  //{
  //  if (drives[i]) total_drives++;
  //}
  //
  //int total_not = nwindows - total_drives;
  //double p = (double) nwindows / (double) total_drives;
  ////p = p*p;
  ////cerr << total_drives << " of " << nwindows << " windows drive expression for penalty of " << p << endl;
  //penalty[&gene] = p;
  
}

void Nuclei::printR2D(Gene& gene, ostream& os)
{
  int window = competition->getWindow();
  int shift  = competition->getShift();
  
  int p = mode->getPrecision();
  int w = p + 7;
  os << gene.getName() << endl;
  competition_data& gcomp = competition_map[&gene];
  
  // print the header
  os << setw(w) << "id";
  for (int i=0; i<gcomp.nwindows; i++)
    os << setw(w) << (i+1)*shift - window/2;
  os << endl;
  
  for (int i=0; i<n; i++)
  {
    os << setw(w) << IDs[i];
    for (int j=0; j<gcomp.nwindows; j++)
      os << setw(w) << setprecision(3) << fixed << gcomp.R_2D[j][i];
    os << endl;
  }
}

void Nuclei::printN2D(Gene& gene, ostream& os)
{
  int window = competition->getWindow();
  int shift  = competition->getShift();
  
  int p = mode->getPrecision();
  int w = p + 7;
  os << gene.getName() << endl;
  competition_data& gcomp = competition_map[&gene];
  
  // print the header
  os << setw(w) << "id";
  for (int i=0; i<gcomp.nwindows; i++)
    os << setw(w) << (i+1)*shift - window/2;
  os << endl;
  
  for (int i=0; i<n; i++)
  {
    os << setw(w) << IDs[i];
    for (int j=0; j<gcomp.nwindows; j++)
      os << setw(w) << setprecision(p) << fixed << gcomp.N_2D[j][i];
    os << endl;
  }
}


  

