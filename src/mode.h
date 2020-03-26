/*********************************************************************************
*                                                                                *
*     mode.h                                                                     *
*                                                                                *
*     Information for how the program is run. If it is something that will       *                     
*     change the output of the program in any way, it should be here so that     *                       
*     a change in output is directly reflected in the input file. Anything       *                     
*     that does not change output (i.e. verbose level) should be handled by      *                      
*     flags, with the possible exception of modes that change how efficiently    *                        
*     the code runs.                                                             *
*                                                                                *
*********************************************************************************/

#ifndef MODE_H
#define MODE_H

#include <map>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/shared_ptr.hpp>

using namespace std;
using boost::property_tree::ptree;

enum annealing_schedules { LAM = 0, EXP = 1};

class Mode
{
private:
  string filename;         // the file name in the command line
  
  double scale_to;         // the value to scale to
  double min_data;         // for functions that divide by data, the minimum to use
  double p_thresh;         // for using p-value thresholds for pwms
  double gc;               // the default gc content to use if unspecified in input
  double penalty_weight;   // the weight given to the penalty function
  double non_specific_k;   // the relative energy of nonspecific binding
  string occupancy_method; // dynamic, combinations, sampling
  string score_function;   // see, chisq, cc, etc
  string scale_data_type;  // what function to use when scaling data
  string seed_string;      // the seed string in the mode
  bool   scale_data;       // do we want to scale the data?
  bool   per_gene;         // report the score per gene
  bool   per_nuc;          // report the score per nuc
  bool   profiling;        // if true, just do initial loop and exit
  bool   self_competition; // whether a TF can compete with itself
  bool   chromatin;        // whether we read in accessibility per gene
  int    verbose;          // how much info to print during running
  int    num_threads;      // if parallel, the number of threads to use
  int    schedule;         // the annealing schhedule to use
  int    precision;        // the precision of printed numbers in output
  unsigned int seed;       // the seed after processing
  int bindingsite_list;   // WSB; get BindingSiteList as K or score (default: PWM)
  
  // promoter competition
  bool   competition; // whether to use this mode
  int    window;      // the window size
  int    shift;       // the shift size
  int    n;           // selectivity coefficient 
  double t;           // threshold for interaction
  
  template< typename T> 
  void readNode(ptree&, string, T*, T);
  void readScaleData(ptree&);
  void readCompetition(ptree&);
  void writeScaleData(ptree&);
  void process_seed_string();
  
public:
  // Constructors
  Mode();
  Mode(string fname, ptree& pt);
  Mode(string fname);
  
  // Getters
  string       getOccupancyMethod()    { return occupancy_method;   }
  string       getScoreFunction()      { return score_function;     }
  string       getScaleDataType()      { return scale_data_type;    }
  double       getScaleTo()            { return scale_to;           }
  double       getMinData()            { return min_data;           }
  double       getPThresh()            { return p_thresh;           }
  double       getGC()                 { return gc;                 }
  double       getPenaltyWeight()      { return penalty_weight;     }
  double       getNonSpecificK()       { return non_specific_k;     }
  bool         getPerGene()            { return per_gene;           }
  bool         getPerNuc()             { return per_nuc;            }
  bool         getProfiling()          { return profiling;          }
  bool         getCompetition()        { return competition;        }
  bool         getSelfCompetition()    { return self_competition;   }
  bool         getScaleData()          { return scale_data;         }
  bool         getChromatin()          { return chromatin;          }
  int          getWindow()             { return window;             }
  int          getShift()              { return shift;              }
  int          getNumThreads()         { return num_threads;        }
  int          getSchedule()           { return schedule;           }
  int          getVerbose()            { return verbose;            }
  int          getPrecision()          { return precision;          }
  int          getN()                  { return n;                  }
  double       getT()                  { return t;                  }
  unsigned int getSeed()               { return seed;               }
  int          getBindingSiteList()    { return bindingsite_list;   }
  
  // Setters
  void setOccupancyMethod(string occupancy_method) { this->occupancy_method = occupancy_method;   }
  void setScoreFunction(string score_function)     { this->score_function   = score_function;     }
  void setScaleDataType(string scale_data_type)    { this->scale_data_type  = scale_data_type;    }
  void setScaleTo(double scale_to)                 { this->scale_to         = scale_to;           }
  void setMinData(double min_data)                 { this->min_data         = min_data;           }
  void setPThresh(double p_thresh)                 { this->p_thresh         = p_thresh;           }
  void setGC(double gc)                            { this->gc               = gc;                 }
  void setPerGene(bool per_gene)                   { this->per_gene         = per_gene;           }
  void setPerNuc(bool per_nuc)                     { this->per_nuc          = per_nuc;            }
  void setProfiling(bool profiling)                { this->profiling        = profiling;          }
  void setCompetition(bool competition)            { this->competition      = competition;        }
  void setSelfCompetition(bool self_competition)   { this->self_competition = self_competition;   }
  void setScaleData(bool scale_data)               { this->scale_data       = scale_data;         }
  void setChromatin(bool chromatin)                { this->chromatin        = chromatin;          }
  void setWindow(int window)                       { this->window           = window;             }
  void setShift(int shift)                         { this->shift            = shift;              }
  void setNumThreads(int num_threads)              { this->num_threads      = num_threads;        }
  void setSchedule(int schedule)                   { this->schedule         = schedule;           }
  void setVerbose(int verbose)                     { this->verbose          = verbose;            }
  void setPrecision(int precision)                 { this->precision        = precision;          }
  void setN(int n)                                 { this->n                = n;                  }
  void setT(int t)                                 { this->t                = t;                  }
  void setSeed(unsigned int seed)                  { this->seed             = seed;               }
  void setBindingSiteList(int bindingsite_list)    { this->bindingsite_list = bindingsite_list;   }
  
  // I/O
  void read(ptree& pt);
  void write(ptree& pt);
  
};


typedef boost::shared_ptr<Mode> mode_ptr;














#endif
