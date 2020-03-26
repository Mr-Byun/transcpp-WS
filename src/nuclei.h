/*********************************************************************************
*                                                                                *
*     nuclei.h                                                                   *
*                                                                                *
*     A nuclei object contains all nuclei which share the same TFs. This means   *
*     they share subgroups and quenching interactions.                           *
*                                                                                *
*     In the original code, many of the calculations on occupancy and quenching  *
*     were multiplications by 0. This class was created to prevent such needless *
*     calculations.                                                              *
*                                                                                *
*********************************************************************************/

#ifndef NUCLEI_H
#define NUCLEI_H 

#include "organism.h"

#include "gene.h"
#include "datatable.h"
#include "promoter.h"
#include "subgroup.h"
#include "quenching.h"
#include "TF.h"
#include "mode.h"
#include "competition.h"
#include "chromatin.h"

/* For the most part, nuclei does not own the private data inside it. It simply
points to the data from it's parent class (Organism). The notable exceptions
are the QuenchingInteractions and Subgroups */

class Organism;

class Nuclei
{
private:
  int            n;
  
  // owned by parent object (organism)
  genes_ptr       genes;
  tfs_ptr         tfs;
  table_ptr       tfdata;
  distances_ptr   distances;
  promoters_ptr   promoters;
  mode_ptr        mode;
  competition_ptr competition;
  chromatin_ptr   chromatin; 
  
  // owned by Nuclei
  bindings_ptr  bindings;
  subgroups_ptr subgroups;
  quenching_ptr quenching;
  modifying_ptr coeffects;
  
  vector<string> IDs;

  map<Gene*, double> penalty;
  
  // if not using promoter competition
  map<Gene*, vector<double> > Ns;
  map<Gene*, vector<double> > Rs;
  
  // if using promoter competition
  bool competition_mode;
  struct competition_data
  {
    int nwindows;
    //X_2D[window_number][nuc_number]
    vector< vector<double> > N_2D;
    vector< vector<double> > T_2D;
    vector< vector<double> > R_2D;
    vector< double> delta_N;
    vector< double> total_N;
  };
  map<Gene*, competition_data> competition_map;
  
  void getRDist(Gene& gene, vector<BindingSite*>& sites, int start_idx, int end_idx, vector<double>& out);
  void calc_pascal(vector<double>);
  void calc_pascal_2D(vector<vector<double> >&, vector<vector<double> >& ret);
  
public:
  // Constructors
  Nuclei();
  Nuclei(Organism* parent, tfs_ptr);
  
  //  Setters
  void setParent(Organism* parent);
  void setTFs(tfs_ptr);
  void setGenes(genes_ptr);
  void setBindings(bindings_ptr);
  void setTFData(table_ptr);
  void setDistances(distances_ptr);
  void setPromoters(promoters_ptr);
  
  void addNuc(string id);
  
  void create(ptree& pt);
  void create();
  void createBindings(ptree& pt);
  void createBindings();
  void createOccupancy();
  void createSubgroups();
  void createCoeffects();
  void createQuenching();
  
  void clear();
  
  //  Getters
  tfs_ptr       getTFs();
  genes_ptr     getGenes();
  bindings_ptr  getBindings();
  table_ptr     getTFData();
  distances_ptr getDistances();
  promoters_ptr getPromoters();
  subgroups_ptr getSubgroups();
  quenching_ptr getQuenching();
  modifying_ptr getCoeffects() { return coeffects; }
  
  vector<double>&          getRate(Gene& gene) {return Rs[&gene];}
  vector<double>&          getN(Gene& gene)    {return Ns[&gene];}
  double&                  getRate(Gene& gene, string& id);
  vector<string>&          getIDs() {return IDs;}
  
  vector< vector<double> >& getR2D(Gene& gene)      { return competition_map[&gene].R_2D; }
  vector< vector<double> >& getN2D(Gene& gene)      { return competition_map[&gene].N_2D; }
  vector< vector<double> >& getT2D(Gene& gene)      { return competition_map[&gene].T_2D; }
  int                       getNWindows(Gene& gene) { return competition_map[&gene].nwindows; }
  int                       getWindow() { return competition->getWindow(); }
  int                       getShift()  { return competition->getShift();  }
  
  int size() {return n;}
  
  // Methods
  bool compareTFs(tfs_ptr); // are the tfs the same as in this object?

  void saveSubgroups()    {subgroups->save();}
  void saveQuenching()    {quenching->save();}
  void saveCoeffects()    {coeffects->save();}
  
  void updateSubgroups()  {subgroups->update();}
  void updateQuenching()  {quenching->update();}
  void updateCoeffects()  {coeffects->update();}
  
  void restoreSubgroups() {subgroups->restore();}
  void restoreQuenching() {quenching->restore();}
  void restoreCoeffects() {coeffects->restore();}
  
  void calcQuenching()    { quenching->calc();}
  void calcCoeffects()    { coeffects->calc();}
  void calcOccupancy()    { subgroups->calc_f();}
  
  // Methods by gene
  void saveSubgroups(Gene& gene)    {subgroups->save(gene);}
  void saveQuenching(Gene& gene)    {quenching->save(gene);}
  void saveCoeffects(Gene& gene)    {coeffects->save(gene);}
  
  void updateSubgroups(Gene& gene)  {subgroups->update(gene);}
  void updateQuenching(Gene& gene)  {quenching->update(gene);}
  void updateCoeffects(Gene& gene)  {coeffects->update(gene);}
  
  void restoreSubgroups(Gene& gene) {subgroups->restore(gene);}
  void restoreQuenching(Gene& gene) {quenching->restore(gene);}
  void restoreCoeffects(Gene& gene) {coeffects->restore(gene);}
  
  void calcQuenching(Gene& gene)    { quenching->calc(gene);}
  void calcCoeffects(Gene& gene)    { coeffects->calc(gene);}
  void calcOccupancy(Gene& gene)    { subgroups->calc_f(gene);}
  
  //void calcN();
  void calcR();
  void calcN(Gene&);
  void calcR(Gene&);
  
  void calcR2(Gene&);
  
  void saveSites(TF& tf)           { bindings->saveSites(tf);          }
  void saveScores(TF& tf)          { bindings->saveScores(tf);         }
  void updateScores(TF& tf)        { bindings->updateScores(tf);       }
  void updateScores()              { bindings->updateScores();         }
  void updateSites(TF& tf)         { bindings->updateSites(tf);        }
  void updateSites()               { bindings->updateSites();          }
  void updateK(TF& tf)             { bindings->updateK(tf);            }
  void updateKandLambda(TF& tf)    { bindings->updateKandLambda(tf);   }
  
  void saveSites(Gene& gene, TF& tf)        { bindings->saveSites(gene, tf);        }
  void saveScores(Gene& gene, TF& tf)       { bindings->saveScores(gene, tf);       }
  void updateScores(Gene& gene, TF& tf)     { bindings->updateScores(gene, tf);     }
  void updateScores(Gene& gene)             { bindings->updateScores(gene);         }
  void updateSites(Gene& gene, TF& tf)      { bindings->updateSites(gene, tf);      }
  void updateSites(Gene& gene)              { bindings->updateSites(gene);          }
  void updateK(Gene& gene, TF& tf)          { bindings->updateK(gene, tf);          }
  void updateKandLambda(Gene& gene, TF& tf) { bindings->updateKandLambda(gene, tf); }
  void updateKandLambda(Gene& gene)         { bindings->updateKandLambda(gene);     }
                                   
  void restoreScores(TF& tf)       { bindings->restoreScores(tf);      }
  void restoreSites(TF& tf)        { bindings->restoreSites(tf);       }
  
  void restoreScores(Gene& gene, TF& tf)       { bindings->restoreScores(gene, tf);    }
  void restoreSites(Gene& gene, TF& tf)        { bindings->restoreSites(gene, tf);     }
                                
  void saveAllOccupancy()          { bindings->saveAllOccupancy();          }
  void restoreAllOccupancy()       { bindings->restoreAllOccupancy();       }
  
  void saveAllOccupancy(Gene& gene)          { bindings->saveAllOccupancy(gene);          }
  void restoreAllOccupancy(Gene& gene)       { bindings->restoreAllOccupancy(gene);       }
  
  void saveModeOccupancy()    { bindings->saveModeOccupancy();    }
  void restoreModeOccupancy() { bindings->restoreModeOccupancy(); }
  
  void saveModeOccupancy(Gene& gene)    { bindings->saveModeOccupancy(gene);    }
  void restoreModeOccupancy(Gene& gene) { bindings->restoreModeOccupancy(gene); }
  
  void saveEffectiveOccupancy()    { bindings->saveEffectiveOccupancy();    }
  void restoreEffectiveOccupancy() { bindings->restoreEffectiveOccupancy(); }
  
  void saveEffectiveOccupancy(Gene& gene)    { bindings->saveEffectiveOccupancy(gene);    }
  void restoreEffectiveOccupancy(Gene& gene) { bindings->restoreEffectiveOccupancy(gene); }
  
  void resizeWindows();
  void resizeWindow(Gene& gene);
  //void updateN() {}
  
  void updateR() {calcR();}
  
  void updateR(Gene& gene) {calcR(gene);}
  double& getPenalty(Gene& gene) { return penalty[&gene]; }
  
  
  // I/O
  void printSubgroups(Gene& g, ostream& os) {subgroups->print(g, os);}
  
  void printOccupancy(Gene& g, ostream& os, bool invert)          {bindings->printTotalOccupancy(g,os, invert);    }
  void printModeOccupancy(Gene& g, ostream& os, bool invert)      {bindings->printModeOccupancy(g,os, invert);     }
  void printEffectiveOccupancy(Gene& g, ostream& os, bool invert) {bindings->printEffectiveOccupancy(g,os, invert);}
  
  void printSites(ostream& os)                     { bindings->printSites(os);}                    
  void printSites(Gene& gene, ostream& os)         { bindings->printSites(gene, os);}             
  void printSites(TF& tf, ostream& os)             { bindings->printSites(tf, os);}            
  void printSites(Gene& gene, TF& tf, ostream& os) { bindings->printSites(gene, tf, os);}
  
  void printScores(Gene& gene, ostream& os)        { bindings->printScores(gene, os); }
  
  void printR2D(Gene& gene, ostream& os);
  void printN2D(Gene& gene, ostream& os);
  
  void write(ptree& output) const                  { bindings->write(output); }
};



typedef boost::shared_ptr<Nuclei>  nuclei_ptr;



#endif
