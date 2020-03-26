/*********************************************************************************
*                                                                                *
*     quenching.h                                                                *
*                                                                                *
*     Contains structure and methods to hold quenching interactions              *
*                                                                                *
*********************************************************************************/

#ifndef QUENCHING_H
#define QUENCHING_H

#include "TF.h"
#include "gene.h"
#include "bindingsite.h"
#include "bindings.h"
#include "distance.h"

#include <boost/shared_ptr.hpp>

using namespace std;

/* Becase we have coactivation and corepression to consider, we need some variables
to decide what repression coefficients and what occupancies to use. The modified or
the regular */

struct QuenchingInteraction
{
  site_ptr actor;
  site_ptr target;
  double   distcoef;
};

typedef map<TF*, map<TF*, vector<QuenchingInteraction> > > gene_quenches;
typedef boost::shared_ptr<gene_quenches>                   gene_quenches_ptr;

class QuenchingInteractions
{
private:
  // quenches[gene][actor][target][index]
  /* if we parallelize over genes the map is not safe with openMP. I have to have
  a map of pointers to maps to prevent the latter from moving around with multiple
  threads */
  map<Gene*, gene_quenches_ptr > quenches;
  map<Gene*, gene_quenches_ptr > saved_quenches;

  genes_ptr     genes;
  tfs_ptr       tfs;
  bindings_ptr  bindings;
  distances_ptr distances;
  
  distance_ptr dist;
  
  //bool hasQuenchingInteractions(Gene&, TF& actor, TF& target);
  void quench_f(vector<double>&,vector<double>&,double);
  
public:
  QuenchingInteractions();
  
  void create(genes_ptr, tfs_ptr, bindings_ptr, distances_ptr);
  
  void set(Gene&, TF& actor, TF& target);
  
  void calc();
  void calc(Gene&);
  void calc(Gene&, TF& actor, TF& target);
  //void calc(int j);
  //void calc(Gene&, TF& actor, TF& target, int j);
  
  void setDistance(Distance& d);
  void setBindings(bindings_ptr b);
  
  void initialize(); // sets mode_occupancy to default;
  void initialize(Gene&); // sets mode_occupancy to default;
  void initialize(Gene&, TF&);
  
  void save();
  void clear();
  void restore();
  void update();
  
  void save(Gene& gene);
  void clear(Gene& gene);
  void restore(Gene& gene);
  void update(Gene& gene);
  
  void printSummary();
};

typedef boost::shared_ptr<QuenchingInteractions> quenching_ptr;

struct ModifyingInteraction
{
  site_ptr     actor;
  site_ptr     target;
  double       distcoef;
  coeffect_ptr coef;
};

typedef map<TF*, map<TF*, vector<ModifyingInteraction> > > gene_mods;
typedef boost::shared_ptr<gene_mods>                       gene_mods_ptr;

class ModifyingInteractions
{
private:
  map<Gene*, gene_mods_ptr > mods;
  map<Gene*, gene_mods_ptr > saved_mods;
  
  genes_ptr     genes;
  tfs_ptr       tfs;
  bindings_ptr  bindings;
  distances_ptr distances;
  
  distance_ptr dist;
  
  void mod_f(vector<double>&,vector<double>&,vector<double>&,double,double);
public:
  ModifyingInteractions();
  
  void create(genes_ptr, tfs_ptr, bindings_ptr, distances_ptr);
  
  void set(Gene&, TF& actor, TF& target, coeffect_ptr);
  
  void calc();
  void calc(Gene&);
  void calc(Gene&, TF& actor, TF& target, coeffect_ptr);

  map<BindingSite*, map<TF*, vector<double> > > getEffq(Gene& gene);
  
  void setDistance(Distance& d);
  void setBindings(bindings_ptr b);

  void initialize(); // sets mode_occupancy to default;
  void initialize(Gene&); // sets mode_occupancy to default;
  void initialize(Gene&, TF&);
  
  void save();
  void clear();
  void restore();
  void update();
  
  void save(Gene& gene);
  void clear(Gene& gene);
  void restore(Gene& gene);
  void update(Gene& gene);
};

typedef boost::shared_ptr<ModifyingInteractions> modifying_ptr;






























#endif


