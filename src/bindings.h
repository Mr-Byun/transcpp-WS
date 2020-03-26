/*********************************************************************************
*                                                                                *
*     bindings.h                                                                 *
*                                                                                *
*     Constains a map of all site on all genes                                   *
*                                                                                *
*********************************************************************************/

#ifndef BINDINGS_H
#define BINDINGS_H

#include "TF.h"
#include "gene.h"
#include "chromatin.h"
#include "bindingsite.h"
#include "datatable.h"

#include <boost/shared_ptr.hpp>

using namespace std;

/* Binding sites were initialized as poitners because it makes it simple to
ensure that subgroups and quenching interactions never point to null. This
may or may not be what we actually want and should be reconsidered in the
future */

typedef map<TF*, site_ptr_vector>          gene_sites_map;
typedef boost::shared_ptr<gene_sites_map>  gene_sites_map_ptr;
typedef map<Gene*, gene_sites_map_ptr>     site_map;
                                       
typedef map<TF*, TFscore>                  gene_scores_map;
typedef boost::shared_ptr<gene_scores_map> gene_scores_map_ptr;
typedef map<Gene*, gene_scores_map_ptr>    scores_map;


class Bindings
{
private:
  int nnuc;
  
  genes_ptr     genes;
  tfs_ptr       tfs;
  table_ptr     tfdata;
  mode_ptr      mode;
  chromatin_ptr chromatin;
  
  vector<string> IDs;
  /* // dont think i need these anymore
  map<int, int> id_2_idx;
  map<int, int> idx_2_id;
  */
  
  //scores[gene][tf] = TFscore
  scores_map scores;
  scores_map saved_scores;
  
  // sites[gene][tf] -> BindingSite
  site_map sites;                    
  site_map saved_sites;
  
  /* it may be useful at some points to access sites according to their
  order on DNA. I have done that here in the Bindings class so that many
  other classes can get this information if necessary */
  map<Gene*, vector<BindingSite*> > ordered_sites_f;
  map<Gene*, vector<BindingSite*> > ordered_sites_r;

  void order_sites(Gene&);
  void add_to_ordered(Gene& gene, TF& tf);
  void eraseTF(Gene& gene, TF& tf);
  void verify_order(Gene& gene, TF& tf);
  
  map<TF*, vector<double> > conc;
  
  //bool hasScores(Gene&, TF&);
  //bool hasSites(Gene&, TF&);
  void createSite(site_ptr_vector& tmp_sites, Gene& gene,
		  TF& tf, int m, int n, double score, double k, char orientation,
		  double kmax, vector<double>& v, int nmodes, double kns);
  void createSite(site_ptr_vector& tmp_sites, Gene& gene, TF& tf,
                  int pos, double bsize, double score, char orientation, 
                  double lambda, double kmax, double maxscore,
                  vector<double>& v, int nmodes, double kns);
  
  void trimOverlaps();
  void trimOverlaps(Gene& gene, TF& tf);
  
public:
  
  // Constructors
  Bindings();
  
  // Getters
  vector<BindingSite*>& getFsites(Gene& gene) { return ordered_sites_f[&gene]; }
  vector<BindingSite*>& getRsites(Gene& gene) { return ordered_sites_r[&gene]; }
  site_ptr_vector& getSites(Gene&, TF&);
  map<TF*, site_ptr_vector>& getSites(Gene&);
  TFscore& getScores(Gene& gene, TF& tf); 
  
  genes_ptr getGenes()  { return genes ; }
  tfs_ptr getTFs()      { return tfs   ; }

  // Setters
  void setMasterBindingIdx();
  void setGenes(genes_ptr g);
  void setTFs(tfs_ptr t)             {tfs       = t; }
  void setTFData(table_ptr c)        {tfdata    = c; }
  void setMode(mode_ptr c)           {mode      = c; }
  void setChromatin(chromatin_ptr c) {chromatin = c; }
  
  void create();
  void clear();
  
  void createScores();
  void createSites();
  void createSites(ptree &pt);
  
  void create(Gene&);
  void clear(Gene&);
  
  void createScores(Gene&);
  void createSites(Gene&);
  
  void createScores(Gene&, TF&);  // score a gene for a tf 
  void createSites(Gene&, TF&);          
  
  void addSite(Gene* g, TF* t, site_ptr b);
  
  // Functions for moving
  void saveScores(TF&);
  void updateScores(TF&);
  void updateScores();
  void restoreScores(TF&);
  
  void saveSites(TF&);
  void updateSites(TF&);
  void updateSites();
  void restoreSites(TF&);
  
  void updateK(TF&);
  void updateKandLambda(TF&);
  
  void saveAllOccupancy();
  void restoreAllOccupancy();
  
  void saveOccupancy();
  void restoreOccupancy();
  
  void saveModeOccupancy();
  void restoreModeOccupancy();
  
  void saveEffectiveOccupancy();
  void restoreEffectiveOccupancy();
  
  // for moving one gene
  void saveScores(Gene&, TF&);
  void updateScores(Gene&, TF&);
  void updateScores(Gene&);
  void restoreScores(Gene&, TF&);
  
  void saveSites(Gene&, TF&);
  void updateSites(Gene&, TF&);
  void updateSites(Gene&);
  void restoreSites(Gene&, TF&);
  
  void updateK(Gene&, TF&);
  void updateKandLambda(Gene&, TF&);
  void updateKandLambda(Gene&);
  
  void saveAllOccupancy(Gene&);
  void restoreAllOccupancy(Gene&);
  
  void saveOccupancy(Gene&);
  void restoreOccupancy(Gene&);
  
  void saveModeOccupancy(Gene&);
  void restoreModeOccupancy(Gene&);
  
  void saveEffectiveOccupancy(Gene&);
  void restoreEffectiveOccupancy(Gene&);
  
  void addNuc(string& nuc_id);
  int  getNnuc() {return nnuc; }
  
  bool isEqual(Bindings& bindings);
  
  void printSites(ostream& os);
  void printSites(Gene& gene, ostream& os);
  void printSites(TF& tf,     ostream& os);
  void printSites(Gene& gene, TF& tf, ostream& os);
  void printSites(Gene& gene, TF& tf, ostream& os, int, int);
  
  void printScores(ostream& os);
  void printScores(Gene& gene, ostream& os);
  
  void printTotalOccupancy(Gene& gene, ostream& os, bool invert);
  void printModeOccupancy(Gene& gene, ostream& os, bool invert);
  void printEffectiveOccupancy(Gene& gene, ostream& os, bool invert);

  void write(ptree& output);
  void write(ptree& pt, BindingSite& b) const;
};




typedef boost::shared_ptr<Bindings>            bindings_ptr;








































#endif


