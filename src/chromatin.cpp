/*********************************************************************************
*                                                                                *
*     chromatin.cpp                                                                *
*                                                                                *
*     Contains the accessibility state for each gene                             *
*                                                                                *
*********************************************************************************/

#include "chromatin.h"

#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <sstream>

#define foreach_ BOOST_FOREACH
#define to_string_ boost::lexical_cast<string>

Chromatin::Chromatin() :
  kacc(double_param_ptr(new Parameter<double>("Kacc","Kacc")))
{}

void Chromatin::read(ptree& parent, genes_ptr g, mode_ptr m)
{
  genes = g;
  mode  = m;
  
  if (mode->getChromatin() == false) return;
  
  ptree& chromatin_node = parent.get_child("Chromatin");
  
  foreach_(ptree::value_type& gene_node, chromatin_node)
  {
    // if this is the parameter kacc
    if (gene_node.first == "Kacc")
      kacc->read( (ptree&) gene_node.second);
    else if (gene_node.first == "Gene")
    {
      ptree& node = gene_node.second;

      string name = node.get<string>("<xmlattr>.name");
      
      Gene& gene = genes->getGene(name);
      int n = gene.length();
      int i = 0;
      accessibility[&gene].resize(n+1);
      vector<double>& gacc = accessibility[&gene];
      
      string acc = node.get<string>("<xmlattr>.data");
      string token;
      
      stringstream s(acc);
      while (getline(s,token,','))
        gacc[i++] = atof(token.c_str());
      
      
      if ( (i-1) != n )
        error("chromatin data (" + to_string_(i-1) + ") must be of same length as gene (" + to_string_(n) + ")");
    }
  }
}

void Chromatin::write(ptree& parent)
{
  if (mode->getChromatin() == false) return;
  
  ptree& chromatin_node = parent.add("Chromatin","");

  ptree & kacc_node = chromatin_node.add("Kacc","");
  kacc->write(kacc_node, mode->getPrecision());
  
  typedef map<Gene*, vector<double> >::iterator it_type;
  
  for(it_type it=accessibility.begin(); it != accessibility.end(); it++)
  {
    Gene& gene = *(it->first);
    ptree & gene_node = chromatin_node.add("Gene","");
    
    gene_node.put("<xmlattr>.name", gene.getName());
    
    stringstream tmp;
    tmp.str("");
    vector<double>& data = it->second;
    int l = data.size();
    //cerr << l << endl;
    //cerr << data[0] << endl;
    tmp << data[0];
    for (int i=1; i<l; i++)
    {
      //cerr << data[i] << endl;
      tmp << "," << data[i];
    }
    
    gene_node.put("<xmlattr>.data", tmp.str());
  }
}

void Chromatin::getParameters(param_ptr_vector& p)
{
  if (mode->getChromatin())
  {
    if (kacc->isAnnealed())
      p.push_back(kacc);
  }
}

void Chromatin::getAllParameters(param_ptr_vector& p)
{
  if (mode->getChromatin())
  {
    p.push_back(kacc);
  }
}

        
      