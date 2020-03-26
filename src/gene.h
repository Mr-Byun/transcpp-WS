/*********************************************************************************
*                                                                                *
*     gene.h                                                                     *
*                                                                                *
*     Contains class definition for genes                                        *
*                                                                                *
*********************************************************************************/

#ifndef GENE_H
#define GENE_H

#include <vector>
#include <list>
#include <string>
#include <deque>
#include <cstdlib>
#include <boost/utility.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <map>

#include "fasta.h"
#include "promoter.h"
#include "twobit.h"
#include "utils.h"
#include "scalefactor.h"
#include "sequence.h"

/* gene class could hold a lot of different information, but for now
it will be pretty simple since we really only care about the sequence, and
position relative to tss */

class Sequence; 

class Gene : boost::noncopyable
{
private:
  /* here i have defined a gene as a sequence with a TSS. It is assumed that
  the input is provided with the TSS pointing in the forward direction */

  string gname;
  string header; // the fasta or 2bit header (chromosome in genome files)
  seq_param_ptr sequence;
  int right_bound;
  int left_bound;
  int tss;
  promoter_ptr     promoter;
  scale_factor_ptr scale;
  bool include;
  double weight; // the weight in the score function
  
  
public:
  // constructors
  Gene(); 
  Gene(string gname, string header, int right_bound, int left_bound, int tss, promoter_ptr p,scale_factor_ptr s, seq_param_ptr seq, double weight);
  
  // getters
  vector<int>&       getSequence();
  const string       getSequenceString() const;
  vector<char>       getSequenceChars() const;
  const string&      getName() const;
  const string&      getChr() const;
  int                getRightBound() const;
  int                getLeftBound() const;
  int                length() const;
  bool               getInclude() { return include; }
  double getRate(double M) { return promoter->getRate(M); } 
  double getWeight() { return weight; }
  scale_factor_ptr   getScale() {return scale;}
  void getParameters(param_ptr_vector& p);
  void getAllParameters(param_ptr_vector& p);
  seq_param_ptr      getSequenceParam() { return sequence; }
  promoter_ptr       getPromoter()      { return promoter; }
  
  // setters
  void setWeight(double weight) { this->weight = weight; }
  void setSequence(vector<int>  s);
  void setSequence(vector<char> s);
  void setSequence(string s);
  void setName(string s);
  void setRightBound(int r);
  void setLeftBound(int l);
  void setInclude(bool include) { this->include = include; }
  
  // methods
  
  // I/O
  void print(ostream& os) const;
  void write(ostream& os) const;
  void write(ptree& pt) const;
};

typedef boost::shared_ptr<Gene> gene_ptr;
typedef vector<gene_ptr> gene_ptr_vector;
typedef boost::shared_ptr<TwoBit> twobit_ptr;
typedef boost::shared_ptr<Fasta>  fasta_ptr;

class GeneContainer
{
private:
  gene_ptr_vector genes;
  scale_factors_ptr scales;
  
  /* this map only serves to remember which genes game from which source
  which leads to more intuitive printing */
  // gene_map[type][sourcename][genes]
  map<string, map<string, gene_ptr_vector> > gene_map; 
  
  /* the source files themselves */
  map<string, twobit_ptr> twobits;
  map<string, fasta_ptr>  fastas;
  promoters_ptr   promoters;
  
  void readTwoBitGenes(ptree& gene_nodes, string& name);
  void readFastaGenes(ptree& gene_nodes, string& name);
  void readLocalGenes(ptree& gene_nodes);
public:
  // Constructor
  GeneContainer();
  GeneContainer(ptree& pt, promoters_ptr, scale_factors_ptr);
  
  // Getters
  Gene& getGene(const string& n);
  Gene& getGene(int index);
  gene_ptr getGeneptr(const string& n);
  gene_ptr getGeneptr(int index);
  int size() const;
  void getParameters(param_ptr_vector& p);
  void getAllParameters(param_ptr_vector& p);
  
  // Setters
  void setPromoters(promoters_ptr p)        { promoters = p;}
  void setScaleFactors(scale_factors_ptr p) { scales = p; }
  void add(gene_ptr);
  
  // I/O
  void read(ptree& pt);
  void write(ostream& os);
  void write(ptree& pt);
};


typedef boost::shared_ptr<GeneContainer> genes_ptr;

#endif


