#include <Rcpp.h>
#include "organism.h"
#include "r_datatable.h"
#include "r_mode.h"
#include "r_defaults.h"
#include "r_gene.h"
#include "r_pwm.h"
#include "r_tf.h"
#include "r_parameter.h"

#ifndef R_ORGANISM_H
#define R_ORGANISM_H

class OrganismPtr
{
private:
  boost::shared_ptr<Organism> organism;
  Defaults par_defaults;
  
public:
  OrganismPtr();
  OrganismPtr(string fname);
  OrganismPtr(string fname, string section);
  
  DataTablePtr    tf_data();
  DataTablePtr    rate_data();
  ModePtr         mode();
  Rcpp::List      get_genes();
  Rcpp::List      tfs();
  Rcpp::List      scaled_rate_data();
  Rcpp::List      bindings();
  Rcpp::List      f();
  Rcpp::List      fa();
  Rcpp::List      F();
  Rcpp::List      rate();
  Rcpp::List      R2D();
  Rcpp::List      N2D();
  Rcpp::List      T2D();
  Rcpp::List      scores();
  Rcpp::List      parameters();
  Rcpp::List      parameter_table();
  double          score()            { return organism->get_score(); }
  Defaults        get_par_defaults() { return par_defaults; }
  CharacterVector get_gene_names();
  CharacterVector get_tf_names();
  CharacterVector get_nuc_names();
  
  void set_genes(Rcpp::List);
  void add_tf(TFPtr);
  
  void            reset_all()   { organism->ResetAll(); }
  void            recalculate() { organism->Recalculate(); }
  
  // dummy functions prevent errors when setting values through organism,
  // i.e "organism$mode$schedule <- 1" attempts to modify mode and schedule,
  // the latter of which is not allowed. If a dummy function is called, it is
  // allowed, but does nothing.
  void listdummy(Rcpp::List x) {}
  void modedummy(ModePtr x)    {}
};
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#endif
