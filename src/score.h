/*********************************************************************************
*                                                                                *
*     score.h                                                                    *
*                                                                                *
*     Presumably we could use any number of functions to score the model.        *
*     We have some set of predictions and outputs, expressed here as pointer     *
*     vectors, as well as some rule for weighting the scores of each             *
*                                                                                *
*********************************************************************************/

#ifndef SCORE_H
#define SCORE_H

#include "organism.h"
#include "gene.h"
#include "TF.h"
#include "bindings.h"
#include "mode.h"
#include "scalefactor.h"

// forward declare organism so that we can get info from it
class Organism;

/*****************  Scoring Functions  ******************************************/

/* Each scoring function takes two vectors of double pointers and outputs a single
double representing the score of that function */

//double sse(vector<double*>& data, vector<double*>& prediction, scale_factor_ptr);

class Score
{
private:
  
  Organism* parent;
  mode_ptr  mode;
  genes_ptr genes;
  
  vector<string> ids;
  
  vector<vector <double*> > data;
  vector<vector <double> >  weighted_data;
  vector<vector <double*> > prediction;
  vector<double*> penalty;

  vector<scale_factor_ptr> scale;
  
  vector<double> scores;
  vector<double> weights;
  
  double scale_to;
  double min_data;
  double divisor;
  
  double score;
  double penalty_weight;
  
  // data scaling functions
  void area();
  void height();
  void none();
  
  void per_gene();
  void per_nuc();
  
  // scoring functions
  typedef void (Score::*SFP)();
  SFP scoreFunc;
  void sse();
  void chisq();
  void pdiff();
  void rms();
  void arkim();
  void cc();
  void rho();
  void sum_slope_squares();
  
  
public:
  
  // Constructors
  Score();
  Score(Organism* parent);
  
  // Getters
  double getScore();
  
  // Setters
  void set(Organism* parent);
  void setWeights();
  
  // Method
  void checkScale(ostream& os);
  
  // I/O
  void print(ostream& os);
  void printMax(ostream& os);
};
  
 
  




#endif
  
  
