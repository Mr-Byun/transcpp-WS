#include <Rcpp.h>
#include <boost/shared_ptr.hpp>
#include "organism.h"

#ifndef R_MODE_H
#define R_MODE_H

class ModePtr
{
private:
  boost::shared_ptr<Mode> mode;
  
public:
  ModePtr();
  ModePtr(string fname);
  ModePtr(string fname, string section);
  ModePtr(mode_ptr mode) { this->mode = mode; } // cannot be exposed!
  
  // getters
  string       get_occupancy_method() { return mode->getOccupancyMethod(); }     
  string       get_score_function()   { return mode->getScoreFunction();   }
  string       get_scale_data_type()  { return mode->getScaleDataType();   }
  double       get_scale_to()         { return mode->getScaleTo();         }
  double       get_min_data()         { return mode->getMinData();         }
  double       get_p_thresh()         { return mode->getPThresh();         }
  double       get_gc()               { return mode->getGC();              }
  bool         get_per_gene()         { return mode->getPerGene();         }
  bool         get_per_nuc()          { return mode->getPerNuc();          }
  bool         get_profiling()        { return mode->getProfiling();       }
  bool         get_competition()      { return mode->getCompetition();     }
  bool         get_self_competition() { return mode->getSelfCompetition(); }
  bool         get_scale_data()       { return mode->getScaleData();       }
  int          get_window()           { return mode->getWindow();          }
  int          get_shift()            { return mode->getShift();           }
  int          get_num_threads()      { return mode->getNumThreads();      }
  int          get_schedule()         { return mode->getSchedule();        }
  int          get_verbose()          { return mode->getVerbose();         }
  int          get_precision()        { return mode->getPrecision();       }
  int          get_n()                { return mode->getN();               }
  unsigned int get_seed()             { return mode->getSeed();            }
  
  // setters
  void set_occupancy_method(string occupancy_method)  { mode->setOccupancyMethod(occupancy_method); }
  void set_score_function(string score_function)      { mode->setScoreFunction(score_function);     }
  void set_scale_data_type(string scale_data_type)    { mode->setScaleDataType(scale_data_type);    }
  void set_scale_to(double scale_to)                  { mode->setScaleTo(scale_to);                 }
  void set_min_data(double min_data)                  { mode->setMinData(min_data);                 }
  void set_p_thresh(double p_thresh)                  { mode->setPThresh(p_thresh);                 }
  void set_gc(double gc)                              { mode->setGC(gc);                            }
  void set_per_gene(bool per_gene)                    { mode->setPerGene(per_gene);                 }
  void set_per_nuc(bool per_nuc)                      { mode->setPerNuc(per_nuc);                   }
  void set_profiling(bool profiling)                  { mode->setProfiling(profiling);              }
  void set_competition(bool competition)              { mode->setCompetition(competition);          }
  void set_self_competition(bool self_competition)    { mode->setSelfCompetition(self_competition); }
  void set_scale_data(bool scale_data)                { mode->setScaleData(scale_data);             }
  void set_window(int window)                         { mode->setWindow(window);                    }
  void set_shift(int shift)                           { mode->setShift(shift);                      }
  void set_num_threads(int num_threads)               { mode->setNumThreads(num_threads);           }
  void set_schedule(int schedule)                     { mode->setSchedule(schedule);                }
  void set_verbose(int verbose)                       { mode->setVerbose(verbose);                  }
  void set_precision(int precision)                   { mode->setPrecision(precision);              }
  void set_n(int n)                                   { mode->setN(n);                              }
  void set_seed(unsigned int seed)                    { mode->setSeed(seed);                        }

};
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#endif
