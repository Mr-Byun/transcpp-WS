/*********************************************************************************
*                                                                                *
*     gene.cpp                                                                   *
*                                                                                *
*     Contains class definition for genes                                        *
*                                                                                *
*********************************************************************************/

#include "gene.h"
#include "utils.h"
#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>

# define foreach_ BOOST_FOREACH


/********************************   Gene    *************************************/


/*    Constructors    */


Gene::Gene() :
  sequence(seq_param_ptr(new Parameter<Sequence>))
{}

Gene::Gene(string name, string head, int left, int right, int t, promoter_ptr p, scale_factor_ptr s, seq_param_ptr seq, double w) 
{
  weight      = w;
  gname       = name;
  header      = head;
  right_bound = right;
  left_bound  = left;
  promoter    = p;
  scale       = s;
  tss         = t;
  sequence    = seq;
  seq->setMove("ResetAll");
  seq->setRestore("ResetAll");
}


/*    Getters   */
  
vector<int>& Gene::getSequence() { return(sequence->getValue().getSequence()); }
  
const string Gene::getSequenceString() const { return(int2string(sequence->getValue().getSequence())); }

vector<char> Gene::getSequenceChars() const 
{
  vector<char> out;
  vector<int>& seq = sequence->getValue().getSequence();
  int length = seq.size();
  out.resize(length);
  
  for (int i=0; i<length; i++)
  {
    switch (seq[i])
    {
    case 0:
      out[i] = 'A';
      break;
    case 1:
      out[i] = 'C';
      break;
    case 2:
      out[i] = 'G';
      break;
    case 3:
      out[i] = 'T';
      break;
    case 4:
      out[i] = 'N';
      break;
    }
  }
  return out;
}

const string & Gene::getName() const { return(gname); }

int Gene::getRightBound() const { return(right_bound); }
    
int Gene::getLeftBound() const {return(left_bound); }

int Gene::length() const { return(sequence->getValue().getLength()); }

void Gene::getParameters(param_ptr_vector& p)
{
  if (sequence->isAnnealed())
    p.push_back(sequence);
}

void Gene::getAllParameters(param_ptr_vector& p)
{
  p.push_back(sequence);
}


/*    Setters   */

void Gene::setSequence(vector<int> s)
{
  sequence->set(s);
  right_bound = 0; // default bound is 0
  left_bound = - length();
}

void Gene::setSequence(vector<char> s)
{
  sequence->set(s);
  right_bound = 0; // default bound is 0
  left_bound = - length();
}

void Gene::setSequence(string s)
{
  sequence->set(s);
  right_bound = 0; // default bound is 0
  left_bound = -length();
}

void Gene::setName(string s) { gname = s; }

void Gene::setRightBound(int r)
{
  right_bound = r;
  left_bound = r-length();
}

void Gene::setLeftBound(int l)
{
  left_bound = l;
  right_bound = l+length();
}


/*    Methods   */


/*    Print   */

void Gene::print(ostream& os) const
{
  os << ">" << gname 
     << "|" << right_bound
     << "|" << length() << endl;
  os << int2string(sequence->getValue().getSequence()) << endl;
}

void Gene::write(ostream& os) const
{
  ptree pt;
  write(pt);
#if BOOST_VERSION / 100 % 1000 < 56
  write_xml_element(os, 
    basic_string<ptree::key_type::value_type>(), 
    pt, 
    -1, 
    boost::property_tree::xml_writer_make_settings<char>(' ', 2));
#else
  write_xml_element(os, 
    basic_string<ptree::key_type::value_type>(), 
    pt, 
    -1, 
    boost::property_tree::xml_writer_make_settings<string>(' ', 2));
#endif
}

void Gene::write(ptree& pt) const
{
  ptree& gene = pt.add("Gene","");
  gene.add("<xmlattr>.name",       gname);
  gene.add("<xmlattr>.header",     header);
  gene.add("<xmlattr>.left_bound", left_bound);
  gene.add("<xmlattr>.right_bound",right_bound);
  gene.add("<xmlattr>.TSS",        tss);
  gene.add("<xmlattr>.promoter",   promoter->getName());
  gene.add("<xmlattr>.weight",     weight);
  gene.add("<xmlattr>.include",    include);
  gene.add("<xmlattr>.sequence", sequence->getValue().getSequenceString());
  
  if (scale->getName() != string("default"))
    gene.add("<xmlattr>.scale",    scale->getName());
}
  

/******************************   GeneContainer    ***********************************/
      
/*    Constructor   */

GeneContainer::GeneContainer() {}

GeneContainer::GeneContainer(ptree& pt, promoters_ptr p, scale_factors_ptr s)
{
  promoters = p;
  scales    = s;
  read(pt);
}


/*    Setters   */

void GeneContainer::add(gene_ptr gene) { genes.push_back(gene); }


/*    Getters   */

Gene& GeneContainer::getGene(int index) {return *genes[index]; }

Gene& GeneContainer::getGene(const string& n) 
{
  int ngenes = genes.size();
  for (int i=0; i<ngenes; i++)
  {
    if (genes[i]->getName() == n)
      return(*genes[i]);
  }
  stringstream err;
  err << "ERROR: getGene() could not find gene with name " << n << endl;
  error(err.str());
  return(*genes[0]);
}

gene_ptr GeneContainer::getGeneptr(int index) {return genes[index]; }

gene_ptr GeneContainer::getGeneptr(const string& n) 
{
  int ngenes = genes.size();
  for (int i=0; i<ngenes; i++)
  {
    if (genes[i]->getName() == n)
      return(genes[i]);
  }
  stringstream err;
  err << "ERROR: getGene() could not find gene with name " << n << endl;
  error(err.str());
  return(genes[0]);
}


int GeneContainer::size() const { return genes.size(); }

void GeneContainer::getParameters(param_ptr_vector& p)
{
  int ngenes = genes.size();
  for (int i=0; i<ngenes; i++)
    genes[i]->getParameters(p);
}

void GeneContainer::getAllParameters(param_ptr_vector& p)
{
  int ngenes = genes.size();
  for (int i=0; i<ngenes; i++)
    genes[i]->getAllParameters(p);
}

/*    I/O    */


void GeneContainer::read(ptree & pt)
{
  ptree& genes_node = pt.get_child("Genes");
  foreach_(ptree::value_type const& source, genes_node)
  {
    if (source.first != "Source") continue;

    string name = source.second.get<string>("<xmlattr>.name", "null");
    string file = source.second.get<string>("<xmlattr>.file", "null");
    string type = source.second.get<string>("<xmlattr>.type", "null");
    if (type == "twobit")
    {
      twobit_ptr tmp(new TwoBit(file));
      twobits[name] = tmp;
      readTwoBitGenes( (ptree&) source.second, name);
    }
    else if (type == "fasta")
    {
      fasta_ptr tmp(new Fasta(file));
      fastas[name] = tmp;
      readFastaGenes( (ptree&) source.second, name);
    }
    else if (type == "local")
    {
       readLocalGenes( (ptree&) source.second );
    }
    else
    {
      stringstream err;
      err << "ERROR: unrecognized sequence source file type " << type << endl;
      error(err.str());
    }
  }
}

void GeneContainer::readTwoBitGenes(ptree& gene_nodes, string& twobit_name)
{
  foreach_(ptree::value_type& gene_node, gene_nodes)
  {
    if (gene_node.first != "Gene") continue;
    
    ptree& node = gene_node.second;
    
    bool include = gene_node.second.get<bool>("<xmlattr>.include", true);
    //if (!include) continue;
    
    int nan = numeric_limits<int>::signaling_NaN();
    
    string gname  = gene_node.second.get<string>("<xmlattr>.name");
    string header = gene_node.second.get<string>("<xmlattr>.header");
    
    int right_bound = gene_node.second.get<int>("<xmlattr>.right_bound", nan); 
    int left_bound  = gene_node.second.get<int>("<xmlattr>.left_bound",  nan);
    int tss         = gene_node.second.get<int>("<xmlattr>.TSS", nan);
    double weight   = gene_node.second.get<double>("<xmlattr>.weight", 1.0);
    
    if (right_bound == nan || left_bound == nan || tss == nan)
    {
      stringstream err;
      err << "ERROR: right_bound, left_bound, and TSS must be set to positive integers when reading from twobit files" << endl;
      error(err.str());
    }
    
    string promoter_name  = gene_node.second.get<string>("<xmlattr>.promoter");
    promoter_ptr promoter = promoters->getPromoter(promoter_name);
    
    string scale_name      = gene_node.second.get<string>("<xmlattr>.scale", "default");
    scale_factor_ptr scale = scales->getScaleFactor(scale_name);
    
    string sequence = gene_node.second.get<string>("<xmlattr>.sequence", "");
    if (sequence == "")
      sequence = twobits[twobit_name]->getSequence(header, left_bound, right_bound);
    
    seq_param_ptr seq(new Parameter<Sequence>());
    seq->setNode(&node);
    seq->set(Sequence(sequence));
    bool anneal = gene_node.second.get<bool>("<xmlattr>.anneal", false);
    seq->setAnnealed(anneal);
    seq->setParamName(gname);

    gene_ptr tmp(new Gene(gname, header, left_bound, right_bound, tss, promoter, scale, seq, weight));
    tmp->setInclude(include);
    genes.push_back(tmp);
    string type;
    if (anneal)
      type = "local";
    else
      type = "twobit";
    
    gene_map[type][twobit_name].push_back(tmp);
  }
}

void GeneContainer::readFastaGenes(ptree& gene_nodes, string& fasta_name)
{
  foreach_(ptree::value_type& gene_node, gene_nodes)
  {
    if (gene_node.first != "Gene") continue;
    
    ptree& node = gene_node.second;
    
    int nan = numeric_limits<int>::signaling_NaN();
    
    bool include = gene_node.second.get<bool>("<xmlattr>.include", true);
    
    string gname  = gene_node.second.get<string>("<xmlattr>.name");
    string header = gene_node.second.get<string>("<xmlattr>.header");
    
    string sequence = gene_node.second.get<string>("<xmlattr>.sequence", "");
    if (sequence == "")
      sequence = fastas[fasta_name]->getSeq(header);

    bool anneal = gene_node.second.get<bool>("<xmlattr>.anneal", false);
    
    seq_param_ptr seq(new Parameter<Sequence>());
    seq->setNode(&node);
    seq->set(Sequence(sequence));
    seq->setAnnealed(anneal);
    seq->setParamName(gname);
    
    int right_bound = gene_node.second.get<int>("<xmlattr>.right_bound", nan); 
    int left_bound  = gene_node.second.get<int>("<xmlattr>.left_bound",  nan);
    // unless otherwise specified, we assume TSS is at position -1
    int tss         = gene_node.second.get<int>("<xmlattr>.TSS",         -1); 
    double weight   = gene_node.second.get<double>("<xmlattr>.weight", 1.0);
    
    // if no left or right bound is specified, make seq immediately left of tss
    if (right_bound == nan && left_bound == nan)
    {
      right_bound = -1; 
      left_bound  = right_bound - sequence.length() + 1;
    } 
    else if (right_bound == nan && left_bound != nan)
    {
      right_bound = left_bound + sequence.length() - 1;
    }
    else // we ignore left_bound if right bound is set 
    {
      left_bound  = right_bound - sequence.length() + 1;
    }
    
    string promoter_name  = gene_node.second.get<string>("<xmlattr>.promoter");
    promoter_ptr promoter = promoters->getPromoter(promoter_name);
    
    string scale_name      = gene_node.second.get<string>("<xmlattr>.scale", "default");
    scale_factor_ptr scale = scales->getScaleFactor(scale_name);
    
    gene_ptr tmp(new Gene(gname, header, left_bound, right_bound, tss, promoter, scale, seq, weight));
    tmp->setInclude(include);
    genes.push_back(tmp);
    string type;
    if (anneal)
      type = "local";
    else
      type = "fasta";
    gene_map[type][fasta_name].push_back(tmp);
  }
}

void GeneContainer::readLocalGenes(ptree& gene_nodes)
{
  foreach_(ptree::value_type& gene_node, gene_nodes)
  {
    if (gene_node.first != "Gene") continue;
    
    ptree& node = gene_node.second;
    
    int nan = numeric_limits<int>::signaling_NaN();
    
    bool include = gene_node.second.get<bool>("<xmlattr>.include", true);
    
    string gname  = gene_node.second.get<string>("<xmlattr>.name");
    string header = gene_node.second.get<string>("<xmlattr>.header");
    
    string sequence = gene_node.second.get<string>("<xmlattr>.sequence");

    seq_param_ptr seq(new Parameter<Sequence>());
    seq->setNode(&node);
    seq->set(Sequence(sequence));
    seq->setAnnealed(gene_node.second.get<bool>("<xmlattr>.anneal", false));
    seq->setParamName(gname);
    
    int right_bound = gene_node.second.get<int>("<xmlattr>.right_bound", nan); 
    int left_bound  = gene_node.second.get<int>("<xmlattr>.left_bound",  nan);
    // unless otherwise specified, we assume TSS is at position -1
    int tss         = gene_node.second.get<int>("<xmlattr>.TSS",         -1); 
    double weight   = gene_node.second.get<double>("<xmlattr>.weight", 1.0);
    
    // if no left or right bound is specified, make seq immediately left of tss
    if (right_bound == nan && left_bound == nan)
    {
      right_bound = -1; 
      left_bound  = right_bound - sequence.length() + 1;
    } 
    else if (right_bound == nan && left_bound != nan)
    {
      right_bound = left_bound + sequence.length() - 1;
    }
    else // we ignore left_bound if right bound is set 
    {
      left_bound  = right_bound - sequence.length() + 1;
    }
    
    string promoter_name  = gene_node.second.get<string>("<xmlattr>.promoter");
    promoter_ptr promoter = promoters->getPromoter(promoter_name);
    
    string scale_name      = gene_node.second.get<string>("<xmlattr>.scale", "default");
    scale_factor_ptr scale = scales->getScaleFactor(scale_name);
    
    gene_ptr tmp(new Gene(gname, header, left_bound, right_bound, tss, promoter, scale, seq, weight));
    tmp->setInclude(include);
    genes.push_back(tmp);
    string type("local");
    gene_map[type]["local"].push_back(tmp);
  }
}
      
void GeneContainer::write(ptree& pt) 
{
  typedef map<string, map<string, gene_ptr_vector> >::iterator it1_type;
  typedef map<string, gene_ptr_vector>::iterator it2_type;
  
  ptree& genes_node = pt.add("Genes","");
  for(it1_type it1=gene_map.begin(); it1 != gene_map.end(); it1++)
  {
    const string& type = it1->first;
    map<string, gene_ptr_vector>& submap = it1->second;
    for(it2_type it2=submap.begin(); it2 !=submap.end(); it2++)
    {
      const string& source_name = it2->first;
      
      ptree& source_node = genes_node.add("Source","");
  
      string file_name;
      if (type == "twobit")
        file_name = twobits[source_name]->getFileName();
      else if (type == "fasta")
        file_name = fastas[source_name]->getFileName();
      else if (type == "local")
        file_name = "local";
      
      source_node.put("<xmlattr>.name", source_name);
      source_node.put("<xmlattr>.file", file_name);
      source_node.put("<xmlattr>.type", type);
      
      int ngenes = gene_map[type][source_name].size();
      for (int i=0; i<ngenes; i++)
        gene_map[type][source_name][i]->write(source_node);
    }
  }   
}

void GeneContainer::write(ostream& os) 
{
  ptree pt;
  write(pt);
#if BOOST_VERSION / 100 % 1000 < 56
  write_xml_element(os, 
    basic_string<ptree::key_type::value_type>(), 
    pt, 
    -1, 
    boost::property_tree::xml_writer_make_settings<char>(' ', 2));
#else
  write_xml_element(os, 
    basic_string<ptree::key_type::value_type>(), 
    pt, 
    -1, 
    boost::property_tree::xml_writer_make_settings<string>(' ', 2));
#endif
}

