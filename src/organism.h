/*********************************************************************************
*                                                                                *
*     organism.h                                                                 *
*                                                                                *
*     An organism is defined here as a collection of nuclei. This is the master  *
*     class which holds most of the data. It creates and array of nuclei         *
*     dynamicaly based on the data present to be scored against.                 *
*     Its main purpose is to be an interface for the model and the annealer.     *
*                                                                                *
*********************************************************************************/

#ifndef ORGANISM_H
#define ORGANISM_H

#include "mode.h"
#include "distance.h"
#include "gene.h"
#include "TF.h"
#include "datatable.h"
#include "score.h"
#include "sequence.h"
#include "nuclei.h"
#include "promoter.h"
#include "parameter.h"
#include "quenching.h"
#include "bindings.h"
#include "scalefactor.h"
#include "coeffects.h"
#include "competition.h"
#include "chromatin.h"

using namespace std;

class   Nuclei;
class   Score;

typedef boost::shared_ptr<Score>  score_ptr;
typedef boost::shared_ptr<Nuclei> nuclei_ptr;

class Organism
{
private:
  mode_ptr          mode;
  distances_ptr     distances;
  tfs_ptr           master_tfs;      // may vary between nuclei
  genes_ptr         master_genes;    // may vary between nuclei
  table_ptr         tfdata;
  table_ptr         ratedata;
  promoters_ptr     promoters;
  scale_factors_ptr scale_factors;
  coeffects_ptr     coeffects;
  score_ptr         score_class;
  coops_ptr         coops;
  competition_ptr   competition;
  chromatin_ptr     chromatin;
  
  param_ptr_vector params; 
  param_ptr_vector all_params;
  
  nuclei_ptr nuclei;
  
  double val;
  double prev;
  
  // we need to translate what is in nuclei objects to an array that corresponds
  // to the embryo itself. When we read in we store the order we read 
  vector<string> ids;
  
  int move_count;
  double thresh; 
  double score_out;
  double previous_score_out;
  
  // the move function that gets called during annealing
  vector<boost::function<void (Organism*)> > moves;
  vector<boost::function<void (Organism*)> > restores;
  
  vector<boost::function<void (Organism*)> > all_moves;
  vector<boost::function<void (Organism*)> > all_restores; 
  
  void setPVectorMoves(vector<boost::function<void (Organism*)> >& mvec, vector<boost::function<void (Organism*)> >& rvec, param_ptr_vector& pvec);
  
  // the move functions
  // void ResetAll(int);  
  void moveScores(TF& tf);
  void movePWM(TF& tf);
  void moveSites(TF& tf);
  void moveLambda(TF& tf);
  void moveKacc();
  void moveKmax(TF& tf);
  void moveCoopD();
  void moveKcoop();
  void moveCoef(double_param_ptr p);
  void moveQuenching();
  void moveQuenchingCoef();
  void moveCoeffect();
  void moveCoeffectEff();
  void movePromoter();
  void moveWindow();
  void null_function();  
  
  // the restore functions
  void restoreScores(TF& tf);
  void restorePWM(TF& tf);
  void restoreSites(TF& tf);
  void restoreLambda(TF& tf);
  void restoreKacc();
  void restoreKmax(TF& tf);
  void restoreCoopD();
  void restoreKcoop();
  void restoreCoef(double_param_ptr p);
  void restoreQuenching();
  void restoreQuenchingCoef();
  void restoreCoeffect(); 
  void restoreCoeffectEff();   

public:
  // Constructors
  Organism();
  Organism(ptree &pt, mode_ptr);
  Organism(string fname, string section);
  
  void initialize(ptree& pt);
  
  // Getters
  table_ptr         getTFData()         {return tfdata;         }
  table_ptr         getRateData()       {return ratedata;       }
  promoters_ptr     getPromoters()      {return promoters;      } 
  distances_ptr     getDistances()      {return distances;      }
  genes_ptr         getGenes()          {return master_genes;   }
  tfs_ptr           getTFs()            {return master_tfs;     }
  mode_ptr          getMode()           {return mode;           }
  scale_factors_ptr getScales()         {return scale_factors;  }
  nuclei_ptr        getNuclei()         {return nuclei;         }
  competition_ptr   getCompetition()    {return competition;    }
  chromatin_ptr     getChromatin()      {return chromatin;      }
  coops_ptr         getCoops()          {return coops;          }
  coeffects_ptr     getCoeffects()      {return coeffects;      }
  vector<string>    getIDs()            {return ids;            }
  double*           getPrediction(Gene&,string&);
  double*           getPenalty(Gene& gene);
  double*           getData(Gene&,string&);
  int               getNNuc()           {return ratedata->getNames("ID").size();}
  int               getNGenes()         {return master_genes->size();}
  double            getTotalScore()     {return get_score();}
  vector<double>&   getN(int gidx);
  vector<double>&   getR(int gidx);
  bindings_ptr      getBindings();
  param_ptr_vector& getParameters()     {return params;} 
  param_ptr_vector& getAllParameters()  {return all_params;}
  int               getNBindings(Gene&);
  
  // Setters
  void populate_nuclei(ptree& pt);
  void populate_nuclei();
  void setMode(mode_ptr mode) {this->mode = mode;}
  void setGenes(genes_ptr genes) {this->master_genes = genes;}
  void addTF(tf_ptr tf) { master_tfs->add(tf); Recalculate(); }
  
  // Methods
  void setMoves();
  void Recalculate(); 
  void ResetAll(); 
  void score();
  void calc_f();
  void scramble();
  void permute(string& table, string& by);
  void checkScale(ostream&);
  void move(int idx);
  void move_all(int idx);
  void restore_all(int idx);
  void adjustThresholds(double percent);

  iparam_ptr getParam(int idx) {return params[idx];}
  string getParamName(int idx);

  template < typename T >
  int  getParamValue(int idx);
  
  // Matlab
  int test_int;
  int test() {test_int++; return test_int;}
  
  // I/O
  void write(string, ptree& pt);
  void printParameters(ostream& os);
  
  void printSites(ostream& os);                     
  void printSites(Gene& gene, ostream& os);         
  void printSites(TF& tf, ostream& os);             
  void printSites(Gene& gene, TF& tf, ostream& os); 
  
  void printScores(Gene& gene, ostream& os);         
  
  void printSubgroups(Gene& gene, ostream& os);
  void printOccupancy(Gene& gene, ostream& os, bool invert);
  void printModeOccupancy(Gene& gene, ostream& os, bool invert);
  void printEffectiveOccupancy(Gene& gene, ostream& os, bool invert);
  
  void printRate(ostream& os, bool invert);
  void printRateData(ostream& os, bool invert);
  void printScore(ostream& os);
  void printMaxScore(ostream& os);
  
  void printR2D(ostream& os);
  void printN2D(ostream& os);
  
  
  // Annealing
  int    getDimension() const; 
  double get_score();
  void   generateMove(int idx, double theta);
  void   restoreMove(int idx);
  void   serialize(void *buf) const;
  void   deserialize(void const *buf);
  int    getStateSize();
};
  

typedef boost::shared_ptr<Organism> organism_ptr;





























#endif
