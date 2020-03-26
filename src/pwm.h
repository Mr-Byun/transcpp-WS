#ifndef PWM_H
#define PWM_H

#include "utils.h"
#include "parameter.h"
#include "mode.h"
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>
#include <vector>

using namespace std;
using boost::property_tree::ptree;

/* there are more ways to express binding preferences than a PWM,
so this is going to have to expand. It could expand to include dinucleotide
preferences, but also more complex rules, like the periodic diculcleotide
preferences of a nulceosome. What they all have in common is that they can
score a sequences, and that score can be converted to an affinity using exp */

/* there are several ways we specify pwm information:

PCM, or Position Count Matrix, contains the counts of 
observed bases at each position

PFM, or Position Frequency Matrix, contains the
probablity of a base at each position

PSSM, the Position Specific Scoring Matrix, contains
the log-likelihood of each base at each position,
taking into account the background model

BEM, the Binding Energy Model, contains what can be
interpreted as the ddG for a substitution of each
base at each position, in the case where lambda=kT
THIS MODE HAS BEEN REMOVED! The conversion of PSSM to BEM
is lossy and I decided there was no real benefit to using
BEM over PSSM. BEMs can still be entered, but they are scored
as PSSMs with a maxscore of 0

*/

/* there are several ways to express binding preferences. We can convert between
some, not others */

enum pwm_types { PCM = 0, PFM = 1, PSSM = 2, BEM = 3, HEIJDEN2012 = 4};

/* struct holding the score of a TF over sequence */
struct TFscore 
{
  //TF* tf;
  vector<double> fscore;
  vector<double> rscore;
  vector<double> mscore;
  double maxscore;
};

class PWM
{
private:
  mode_ptr mode;
  string source;
  int    input_type; // the type of binding preference initially specified
  double gc;         // gc content used as background model
  string tfname;     // the tfname so we know what to move in parameters
   
  // pwm specific parameters
  bool              is_pwm;          // is this even a pwm
  pwm_param_ptr     mat;             // the actual matrix
  vector<int>       consensus;       // the consensus sequence
  vector<double>    position_counts;
  vector<double>    s2p;             // convert scores to pvalues
  double            pseudo;          // pseudo count if PCM
                   
  // nucleosome binding parameters
  bool                    is_periodic; // uses periodic dinuc binding from van der Heijden 2012
  double                  plength;     // the length of sequence to scan over (use 147 or 146 for octamer, 74 for tetramer)
  double_param_ptr        period;      // the period of dinuc (10-11)
  double_param_ptr        beta;        // the strengh of preference for dinucs
  
  double maxscore; // the maximum score

  // private methods
  void subscore(const vector<int> & s, double * out);
  double score_dyad(int first, int second, double position);
  
public:
  // constructors
  PWM();
  PWM(mode_ptr mode);
  PWM(vector<vector<double> >& t, int type, mode_ptr mode);
  
  // setters
  void setSource(string source) { this->source = source; }
  void setGC(double gc);         
  void setPseudo(double pseudo);
  void setTFName(string tfname) { this->tfname = tfname; }
  void setMode(mode_ptr mode) { this->mode = mode; }
  
  void setPWM(vector<vector<double> >& t, int type); // use default gc=0.25, pseudo=1
  void setPWM(vector<vector<double> >& t, int type, double gc, double pseudo);

  // getters
  const string& getSource()    { return source; }
  double getGC()        { return gc; }
  double getPseudo()    { return pseudo; }
  double getMaxScore()  { return maxscore; }
  void setMaxScore(double mscore)
                        { maxscore = mscore; } // WSB
  int    length()       { return mat->getValue().size(); }
  int    getInputType() { return input_type; }
  vector<int>& getConsensus() { return consensus; }
  vector<vector<double> >& getPWM(); // efficient, return reference to pwm
  vector<vector<double> >  getPWM(int type); // return a copy, for conveneince, not inner loop
  pwm_param_ptr getPWMParam() { return mat; }
  void getParameters(param_ptr_vector&);
  void getAllParameters(param_ptr_vector&);

  // forward converstions
  void PCM2PFM();
  void PFM2PSSM();

  // reverse conversions (dont touch the actual matrix since PSSMs are always used internally)
  void PFM2PCM(vector<vector<double> >& t);
  void PSSM2PFM(vector<vector<double> >& t);

  // methods
  void   setNscore();
  void   calc_max_score();
  double pval2score(double pval);  // returns the threshold that would yeild a given p-value
  double score2pval(double score); // returns the pvalue of a given score
  void   score(const vector<int>& s, TFscore &t);
  
  //size_t getSize();
  //void   serialize(void *buf) const;
  //void   deserialize(void const *buf);
  
  //void print(ostream& os, int precision);
  void read(ptree& pt);
  void write(ptree& pt);
};

#endif

