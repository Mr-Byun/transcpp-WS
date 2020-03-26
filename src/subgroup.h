/*********************************************************************************
*                                                                                *
*     subgroup.h                                                                 *
*                                                                                *
*     Contains the subgroup classes and methods.                                 *
*                                                                                *
*********************************************************************************/

#ifndef SUBGROUP_H
#define SUBGROUP_H

#include "bindingsite.h"
#include "bindings.h"
#include "TF.h"
#include "gene.h"
#include "mode.h"
#include <list>
#include <map>
#include <bitset>


/* for every site we need to know some information for calculating either
the full partition function, or parts of it */
/*  Note: the partition index at 0 is the null bindings state, so the partition
index is the site index + 1 */

struct Partition 
{
  bool coops; // flag this site as to whether it cooperates at all
  int last;   // the last partition index that does not compete with the current
  vector<int>      coop_site; // the site index that coops with this
  vector<int>      coop_past; // the last parition index that the coop state doesnt compete with;
  vector<double>   dist_coef; // the distance coefficient of this cooperative interaction
  vector<coop_ptr> coop;
  
#ifdef LARGENUMS
  vector<long double> Z;   // the partition function of sites through partition index
  vector<long double> Zc;  // the partial partition function, cooperating
  vector<long double> Znc; // the partial partition function, non-cooperating 
#elif defined VERYLARGENUMS
  vector<double> log_Z;      // the partition function of sites through partition index
  vector<double> log_Zc;     // the partial partition function, cooperating
  vector<double> log_Znc;    // the partial partition function, non-cooperating 
#else
  vector<double> Z;   // the partition function of sites through partition index
  vector<double> Zc;  // the partial partition function, cooperating
  vector<double> Znc; // the partial partition function, non-cooperating 
#endif
  
};
  


/* A subgroup is a collection of pointers to binding sites. */

class Subgroup
{
private:
  vector<BindingSite*> sites_f; // the pointers to sites in the subgroup
  vector<BindingSite*> sites_r; // the pointers to sites in the subgroup, in reverse order
  vector<int> f2r;              // maps forward index onto reverse index
  
  int right_bound;            // the right most position in the set
  int left_bound;             // the left most position in the set
  bindings_ptr  bindings;     // pointer to master bindings

  vector<Partition> ZF; // forward partition function
  vector<Partition> ZR; // reverse partition function
  
  void pre_process_pair(vector<Partition>&, int, BindingSite*, TF*, char, int, BindingSite*, TF*, char);
  void iterate_partition(vector<Partition>& Z, vector<BindingSite*>& sites, int site_index);

public:
  // constructors
  Subgroup();
  Subgroup(BindingSite*, bindings_ptr);
  
  // setters
  void addSite(BindingSite*);
  void addSubgroup(Subgroup&);
  
  // getters
  int                   size() const;
  vector<BindingSite*>& getSites();
  int                   getRightBound() const;
  int                   getLeftBound() const;
  
  // methods
  bool overlaps(BindingSite&);
  bool overlaps(BindingSite*, BindingSite*);
  void sort();
  bool checkCoop(BindingSite& site);
  
  void occupancy();

  void pre_process();

  // print
  void print(ostream& os);
};

  
typedef boost::shared_ptr<list<Subgroup> > gene_groups_ptr;
class Subgroups
{
private:
  //maps gene to subgroup
  map<Gene*, gene_groups_ptr > groups;
  map<Gene*, gene_groups_ptr > saved_groups;
  
  genes_ptr     genes;
  tfs_ptr       tfs;
  bindings_ptr  bindings;
  mode_ptr      mode;
  
  void addSites(Gene&);
  void addSite(list<Subgroup>&, site_ptr);
  
public:
  Subgroups();
  
  void clear();
  void save();
  void update();
  void restore();
  
  void clear(Gene&);
  void save(Gene&);
  void update(Gene&);
  void restore(Gene&);
  
  void create(genes_ptr, tfs_ptr, bindings_ptr, mode_ptr); 

  void calc_f();
  void calc_f(Gene&);
  //void calc_f(int nuc_idx);
  
  void print(ostream& os);
  void print(Gene&, ostream& os);
  
};
  
  
typedef boost::shared_ptr<Subgroups> subgroups_ptr;
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#endif
