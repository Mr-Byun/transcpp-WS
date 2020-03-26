/*********************************************************************************
*                                                                                *
*     TF.cpp                                                                     *
*                                                                                *
*     Contains class definition for Transcription Factos                         *
*                                                                                *
*********************************************************************************/

# include "TF.h"
# include "utils.h"

# include <cmath>
# include <iomanip>
# include <sstream>
# include <boost/property_tree/xml_parser.hpp>
# include <boost/property_tree/ptree.hpp>
# include <boost/foreach.hpp>
# include <boost/version.hpp>



using namespace std;
using boost::property_tree::ptree;

# define foreach_ BOOST_FOREACH


/*********************************   TF    **************************************/

/*      Constructors    */

TF::TF():
  kmax(double_param_ptr(new Parameter<double>("Kmax","Kmax"))),
  threshold(double_param_ptr(new Parameter<double>("Threshold","Sites"))),
  lambda(double_param_ptr(new Parameter<double>("Lambda","Lambda")))
{
  // initialize all the member variables
  index  = 0;
  bsize  = 0;
  offset = 0;
  Kns    = 0;
  double_param_ptr coef(new Parameter<double>());
  coef->setParamName("coef");
  coef->set(0.0);
  coef->setTFName(tfname);
  coef->setAnnealed(false);
  coef->setMove("Coef");
  coef->setRestore("Coef");
  coef->setType("double");
  coef->setLimits(-1.0,1.0);
  coefs.push_back(coef);
}

TF::TF(ptree& pt, mode_ptr m):
  kmax(double_param_ptr(new Parameter<double>("Kmax","Kmax"))),
  threshold(double_param_ptr(new Parameter<double>("Threshold","Sites"))),
  lambda(double_param_ptr(new Parameter<double>("Lambda","Lambda")))
{
  mode = m;
  set(pt);
}

void TF::set(ptree& pt)
{
  if (!pt.count("PWM")) error("TF::TF(ptree& pt) TF node malformed");

  tfname     = pt.get<string>("<xmlattr>.name");
  bsize      = pt.get<int>("<xmlattr>.bsize");
  
  Kns = mode->getNonSpecificK();
  
  // read in coefficients
  ptree& coefs_node = pt.get_child("Coefficients");
  foreach_(ptree::value_type const& coef_node, coefs_node)
  {
    if (coef_node.first != "coef") continue;
    
    double_param_ptr coef(new Parameter<double>());
    coef->read( (ptree&) coef_node.second);
    coef->setParamName("coef");
    coef->setTFName(tfname);
    coefs.push_back(coef);
  }
  
  ptree& kmax_node = pt.get_child("kmax");
  kmax->read(kmax_node);
  kmax->setTFName(tfname);
  
  ptree& threshold_node = pt.get_child("threshold");
  if(mode->getBindingSiteList())
    {
      threshold->setAnnealed(false);
      threshold->setTFName(tfname);
    }
  // WSB; threshold is disabled during BindingSiteList Mode due to some difficulties in annealing
  else
    {
      threshold->read(threshold_node);
      threshold->setTFName(tfname);
    }
  
  ptree& lambda_node = pt.get_child("lambda");
  if(mode->getBindingSiteList() == 1)
    {
      lambda->setAnnealed(false);
      lambda->setTFName(tfname);
    }
  else
    {
      lambda->read(lambda_node);
      lambda->setTFName(tfname);
    }
  
  readPWM(pt);
  
  double pThresh = mode->getPThresh();
  if (pThresh)
  {
    threshold->set(energy.pval2score(pThresh));
    threshold->setAnnealed(false);
    //cerr << pThresh << " p-value cuttoff for " << tfname << " is " << threshold->getValue() << endl; 
  }
    
  
}

void TF::readPWM(ptree& tf_node)
{
  energy.setMode(mode);
  energy.setTFName(tfname);
  energy.read(tf_node);
  
}
    



/*    Setters   */

void TF::setName(string n) { tfname = n; }

void TF::setPWM(vector<vector<double> >& t, int type) // use default gc=0.25, pseudo=1
{
  energy.setPWM(t, type);
}

void TF::setPWM(vector<vector<double> >& t, int type, double gc, double pseudo)
{
  energy.setPWM(t, type, gc, pseudo);
}
 
void TF::setKmax(double k) { kmax->set(k); }

void TF::setLambda(double l) { lambda->set(l); }

void TF::setThreshold(double t) {threshold->set(t); }

void TF::setBindingSize(int b) { bsize = b; }
  
void TF::setCoefs(vector<double> c)
{
  int input_size = c.size();
  
  for (int i=0; i<input_size; i++)
  {
    if (i < (int) coefs.size())
      coefs[i]->set(c[i]);
    else
    {
      double_param_ptr coef(new Parameter<double>());
      coef->setParamName("coef");
      coef->set(c[i]);
      coef->setTFName(tfname);
      coef->setAnnealed(false);
      coef->setMove("Coef");
      coef->setRestore("Coef");
      coef->setType("double");
      coef->setLimits(-1.0,1.0);
      coefs.push_back(coef);
    }
  }
}
    
/*    Getters   */

const string& TF::getName() const { return tfname; }

double TF::getThreshold() { return threshold->getValue(); }

double TF::getKmax() { return kmax->getValue(); }

double& TF::getCoef() { return coefs[0]->getValue(); }

double TF::getModifiedCoef() 
{
  if (coefs.size() > 1)
    return coefs[0]->getValue();
  else
    return 0;
}

bool TF::neverActivates()
{
  int ncoefs = coefs.size();
  for (int i=0; i<ncoefs; i++)
  {
    if (coefs[i]->getValue() > 0)
      return false;
  }
  return true;
}

bool TF::neverQuenches()
{
  int ncoefs = coefs.size();
  for (int i=0; i<ncoefs; i++)
  {
    if (coefs[i]->getValue() < 0)
      return false;
  }
  return true;
}

vector<double> TF::getCoefs()
{
  vector<double> out;
  int ncoefs = coefs.size();
  out.resize(ncoefs);
  for (int i=0; i<ncoefs; i++)
    out[i] = coefs[i]->getValue();
  return out;
}

double TF::getLambda() { return lambda->getValue(); }

double TF::getMaxScore()  
{ 
  return energy.getMaxScore(); 
}

int TF::getBindingSize() const { return bsize; }

void TF::getParameters(param_ptr_vector& p)
{
  if (kmax->isAnnealed())
    p.push_back(kmax);
  if (threshold->isAnnealed())
    p.push_back(threshold);
  if (lambda->isAnnealed())
    p.push_back(lambda);
  
  int ncoefs = coefs.size();
  for (int i=0; i<ncoefs; i++)
  {
    if (coefs[i]->isAnnealed())
      p.push_back(coefs[i]);
  }
  
  energy.getParameters(p);

}

void TF::getAllParameters(param_ptr_vector& p)
{
  p.push_back(kmax);
  p.push_back(threshold);
  p.push_back(lambda);
  
  int ncoefs = coefs.size();
  for (int i=0; i<ncoefs; i++)
    p.push_back(coefs[i]);

  energy.getAllParameters(p);
  
}

/*bool TF::checkCoops(TF* t)
{
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
  {
    if (coops[i].first == t)
      return true;
  }
  return false;
}*/

// check if cooporates with orientation
bool TF::checkCoops(TF* t, char o1, char o2) 
{
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
  {
    if (coops[i].first == t)
    {
      coop_ptr coop = coops[i].second;
      if (o1 == 'F')
      {
        if      (o2 == 'F' && coop->getHT())
          return true;
        else if (o2 == 'R' && coop->getHH())
          return true;
      } 
      else //(o1 == 'R')
      {
        if      (o2 == 'F' && coop->getTT())
          return true;
        else if (o2 == 'R' && coop->getHT())
          return true;
      }
    }
  }
  return false;
}

coop_ptr TF::getCoop(TF* t)
{
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
  {
    if (coops[i].first == t)
      return coops[i].second;
  }
  coop_ptr void_coop;
  return void_coop;
}
  

/*    Methods   */

void TF::score(const vector<int>& s, TFscore &t)
{
  //PWM& pwm = energy->getValue();
  energy.score(s,t);
}

    
TFscore TF::score(const vector<int>& s)
{
  TFscore t;
  score(s,t);
  return(t);
}

TFscore TF::score(const string & s)
{
  
  vector<int> out(s.size());
  out = string2int(s);
  return(score(out));
}
  

/*    print   */


void TF::write(ostream& os) 
{
  ptree pt; // boost property tree
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


void TF::write(ptree& tfsnode) 
{
  int p = mode->getPrecision();
  stringstream tmp;
  tmp << setprecision(p);
  
  ptree & tfnode = tfsnode.add("TF","");
  tfnode.put("<xmlattr>.name",      tfname);
  
  ptree & kmax_node      = tfnode.add("kmax          ","");
  kmax->write(kmax_node, mode->getPrecision());
  
  ptree & threshold_node  = tfnode.add("threshold     ","");
  threshold->write(threshold_node, mode->getPrecision());

  ptree & lambda_node     = tfnode.add("lambda        ","");
  lambda->write(lambda_node, mode->getPrecision());
  
  ptree& coefs_node = tfnode.add("Coefficients", "");
  int ncoefs = coefs.size();
  for (int i=0; i<ncoefs; i++)
  {
    ptree& coef_node = coefs_node.add("coef","");
    coefs[i]->write(coef_node, mode->getPrecision());
  }
  
  tfnode.put("<xmlattr>.bsize",   bsize);
  tfnode.put("<xmlattr>.include", "true");
  
  energy.write(tfnode);
}
  
    

/*********************************   TFContainer    *********************************/


/*    Constructor   */

TFContainer::TFContainer() {}

TFContainer::TFContainer(ptree& pt, mode_ptr m) 
{
  mode = m;
  add(pt, mode);
}
    

/*    Getters   */

TF& TFContainer::getTF(const string& n)
{
  int i;
  int ntfs = tfs.size();
  
  for (i=0; i<ntfs; i++)
  {
    if(tfs[i]->getName() == n) // they are equal
      return(*tfs[i]);                // return this tf
  }
  // if we are here TF was not found
  stringstream err;
  err << "ERROR: getTF() TF with name " << n << " not found" << endl;
  error(err.str());
  return *tfs[0]; // you will never get here!
}

int TFContainer::size() const { return(tfs.size()); }

TF& TFContainer::getTF(int index) { return *tfs[index]; }

tf_ptr TFContainer::getTFptr(int index) { return tfs[index]; }

tf_ptr TFContainer::getTFptr(const string& n)
{
  int i;
  int ntfs = tfs.size();
  
  for (i=0; i<ntfs; i++)
  {
    if(tfs[i]->getName() == n) // they are equal
      return(tfs[i]);                // return this tf
  }
  // if we are here TF was not found
  stringstream err;
  err << "ERROF: getTFptr() TF with name " << n << " not found" << endl;
  error(err.str());
  return tfs[0]; // you will never get here!
}
  
void TFContainer::getParameters(param_ptr_vector& p)
{
  int ntfs = tfs.size();
  for (int i=0; i<ntfs; i++)
    tfs[i]->getParameters(p);
}

void TFContainer::getAllParameters(param_ptr_vector& p)
{
  int ntfs = tfs.size();
  for (int i=0; i<ntfs; i++)
    tfs[i]->getAllParameters(p);
}
  
  
/*  Setters   */

void TFContainer::setCoops(coops_ptr c)
{
  int ntfs = tfs.size();
  for (int i=0; i<ntfs; i++)
  {
    coop_pairs cur_coops = c->getCoops(tfs[i]->getName());
    // string comparisons are slow, so we want the final object to be TF pointers
    vector< pair<TF*,coop_ptr> > tmp_pairs;
    int ncur_coops = cur_coops.size();
    for (int j=0; j<ncur_coops; j++)
    {
      TF& tf_ref = getTF(cur_coops[j].first);
      TF* tf = &tf_ref;
      tmp_pairs.push_back(pair<TF*,coop_ptr>(tf, cur_coops[j].second));
    }
    tfs[i]->setCoops(tmp_pairs);
  }
}

void TFContainer::setCoeffects(coeffects_ptr c)
{
  int ntfs = tfs.size();
  for (int i=0; i<ntfs; i++)
  {
    coeffect_pairs cur_coefs = c->getTargets(tfs[i]->getName());
    // string comparisons are slow, so we want the final object to be TF pointers
    vector< pair<TF*, coeffect_ptr> > tmp_pairs;
    int ncur_coefs = cur_coefs.size();
    for (int j=0; j<ncur_coefs; j++)
    {
      TF& tf_ref = getTF(cur_coefs[j].first);
      TF* tf     = &tf_ref;
      tmp_pairs.push_back(pair<TF*,coeffect_ptr>(tf,cur_coefs[j].second));
    }
    tfs[i]->setCoeffects(tmp_pairs);
  }
} 
  
  
// can handle getting a TF node or TFs node for single or multiple add
void TFContainer::add(ptree& pt, mode_ptr m)
{
  mode = m;
  if (pt.count("TFs")) // we are at a parent of TFs
  {
    ptree& tfs_node = pt.get_child("TFs");
    foreach_(ptree::value_type const& tf, tfs_node)
    {
      if (tf.first == "TF")
      {
        ptree& tf_node = (ptree&) tf.second;
        if (tf_node.get<bool>("<xmlattr>.include"))
        {
          tf_ptr t(new TF(tf_node, mode) );
          tfs.push_back(t);
          int ntfs = tfs.size();
          tfs[ntfs-1]->setIndex(ntfs);
        }
      }
    }
  }
  else if (pt.count("TF")) // we are at the TFs node
  {
    foreach_(ptree::value_type const& tf, pt)
    {

      if (tf.first == "TF")
      {
        ptree& tf_node = (ptree&) tf.second;
        if (tf_node.get<bool>("<xmlattr>.include"))
        {
          tf_ptr t(new TF(tf_node, mode) );
          tfs.push_back(t);
          int ntfs = tfs.size();
          tfs[ntfs-1]->setIndex(ntfs);
        }
      }
    }
  } 
  else // this is a TF node
  {
    if (pt.get<bool>("<xmlattr>.include"))
    {
      tf_ptr t(new TF(pt, mode) );
      tfs.push_back(t);
      int ntfs = tfs.size();
      tfs[ntfs-1]->setIndex(ntfs);
    }
  }
} 

void TFContainer::add(tf_ptr t) {tfs.push_back(t);}

/*    Output    */
void TFContainer::write(ostream& os) const
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

void TFContainer::write(ptree& pt) const
{
  int i;
  int ntfs = tfs.size();
  ptree & tfsnode = pt.add("TFs","");
  for (i=0; i<ntfs; i++)
    tfs[i]->write(tfsnode);
  
}
   

