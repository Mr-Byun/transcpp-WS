
#include "pwm.h"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <limits.h>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>
#define to_ boost::lexical_cast
#define to_string_ boost::lexical_cast<string>
#define foreach_ BOOST_FOREACH

PWM::PWM():
  mat(pwm_param_ptr(new Parameter<vector<vector<double> > >("PWM","ResetAll"))),
  period(double_param_ptr(new Parameter<double>("Period","ResetAll"))),
  beta(double_param_ptr(new Parameter<double>("Beta","ResetAll")))
{
  // default parameters
  pseudo     = 1.0;
  gc         = 0.5;
  maxscore   = 0;
  input_type = -1;
  is_pwm     = true;
  source     = string("");
}

PWM::PWM(mode_ptr mode):
  mat(pwm_param_ptr(new Parameter<vector<vector<double> > >("PWM","ResetAll"))),
  period(double_param_ptr(new Parameter<double>("Period","ResetAll"))),
  beta(double_param_ptr(new Parameter<double>("Beta","ResetAll")))
{
  // default parameters
  pseudo     = 1.0;
  gc         = 0.5;
  maxscore   = 0;
  input_type = -1;
  this->mode = mode;
  is_pwm     = true;
  source     = string("");
}

PWM::PWM(vector<vector<double> >& t, int type, mode_ptr mode):
  mat(pwm_param_ptr(new Parameter<vector<vector<double> > >("PWM","ResetAll"))),
  period(double_param_ptr(new Parameter<double>("Period","ResetAll"))),
  beta(double_param_ptr(new Parameter<double>("Beta","ResetAll")))
{
  // default parameters
  this->mode = mode;
  source     = string("");
  pseudo     = 1.0;
  gc         = 0.5;
  maxscore   = 0;
  input_type = type;
  is_pwm     = true;

  setPWM(t, type);
}

void PWM::read(ptree& pt)
{
  ptree& pwm_node = pt.get_child("PWM");
  
  source = pwm_node.get<string>("<xmlattr>.source","");
  gc     = pwm_node.get<double>("<xmlattr>.gc", mode->getGC());
  string type = pwm_node.get<string>("<xmlattr>.type");
  
  if (type == "Heijden2012")
  {
    input_type = HEIJDEN2012;
    is_pwm      = false;
    is_periodic = true;
    plength      = pwm_node.get<double>("<xmlattr>.length", 147);
    ptree& beta_node   = pwm_node.get_child("Beta");
    ptree& period_node = pwm_node.get_child("Period");
    
    beta->read(beta_node);
    beta->setTFName(tfname);
    
    period->read(period_node);
    period->setTFName(tfname);
    calc_max_score();
  }
  else
  {
    is_pwm      = true;
    is_periodic = false;
    vector<vector<double> > matrix;
    
    mat->setAnnealed(pwm_node.get<bool>("<xmlattr>.anneal", false));
    mat->setParamName(tfname + "_pwm");
    mat->setNode(&pwm_node);
    mat->setTFName(tfname);
      
    string         token;
    vector<double> line;
    foreach_(ptree::value_type const& v, pwm_node)
    {
      if (v.first == "position")
      {
        line.clear();
        stringstream s(v.second.data());
        while( getline(s, token, ';'))
          line.push_back(atof(token.c_str()));
        line.push_back(0);
        matrix.push_back(line);
      }
    }
    if (type == "PCM")
      setPWM(matrix, PCM);
    else if (type == "PFM")
      setPWM(matrix, PFM);
    else if (type == "PSSM")
      setPWM(matrix, PSSM);
    else if (type == "BEM")
      setPWM(matrix, PSSM);
    else
      error("ERROR: readPWM() did not recognize pwm of type " + type);
    setNscore();
  }
}

void PWM::write(ptree& tfnode)  
{
  ptree & pwmnode = tfnode.add("PWM","");
  if (mat->isAnnealed())
    input_type = PSSM;
  
  switch (input_type)
  {
  case PCM:
    pwmnode.put("<xmlattr>.type","PCM");
    pwmnode.put("<xmlattr>.pseudo", getPseudo());
    break;
  case PFM:
    pwmnode.put("<xmlattr>.type","PFM");
    break;
  case PSSM:
    pwmnode.put("<xmlattr>.type","PSSM");
    break;
  case BEM:
    pwmnode.put("<xmlattr>.type","PSSM");
    break;
  case HEIJDEN2012:
    pwmnode.put("<xmlattr>.type","Heijden2012");
    break;
  default:
    error("TF::write() unrecognized pwm type");
    break;
  }
    
  if (mat->isAnnealed())
    pwmnode.put("<xmlattr>.anneal", "true");
  if (getSource() != string(""))
    pwmnode.put("<xmlattr>.source", getSource());
  if (getGC() != mode->getGC())
    pwmnode.put("<xmlattr>.gc", getGC());
  
  if (input_type != HEIJDEN2012)
  {
    int p = mode->getPrecision();
    stringstream tmp;
    tmp << setprecision(p);
    int w = p + 7;
    
    tmp.str("");
    tmp << setw(w+1) << "A"
        << setw(w+1) << "C"
        << setw(w+1) << "G"
        << setw(w+1) << "T";
    pwmnode.add("base   ", tmp.str());
    vector< vector<double> > matrix = getPWM(input_type);
    int pwmlen=matrix.size();
    
    for (int i=0; i<pwmlen; i++)
    {
      tmp.str("");
      if (input_type == PCM)
      {
        int zeros = (int) pow(10.0,p);
        for (int j=0; j<4; j++)
        {
          if (matrix[i][j] < 0)
            matrix[i][j] = 0.0;
          matrix[i][j] = roundf(matrix[i][j] * zeros) / zeros;
        }
      }
      tmp << setw(w) << matrix[i][0] << ";"
          << setw(w) << matrix[i][1] << ";"
          << setw(w) << matrix[i][2] << ";"
          << setw(w) << matrix[i][3];
      pwmnode.add("position",tmp.str());
    }
  }
  else
  {
    ptree& beta_node = pwmnode.add("Beta   ", "");
    beta->write(beta_node, mode->getPrecision());
    ptree& period_node = pwmnode.add("Period", "");
    period->write(period_node, mode->getPrecision());
    pwmnode.put("<xmlattr>.length", plength);
  }
  
}

void PWM::setPWM(vector<vector<double> >& t, int type, double g, double p)
{
  pseudo     = p;
  gc         = g;
  input_type = type;
  maxscore   = 0;
  
  setPWM(t, type);
}

void PWM::setPWM(vector<vector<double> >& e, int type)
{
  vector<vector<double> >& matrix = mat->getValue();
  matrix = e;
  input_type = type;
  maxscore   = 0;
  
  // we need to verify the integrity of this pwm
  int pwmlen=matrix.size();
  if (pwmlen <= 0)
    error("ERROR: setPWM() with matrix size less than 1");
  for (int i=0; i<pwmlen; i++)
  { 
    int n = matrix[i].size();
    if (n == 4)
      matrix[i].push_back(0);
    else if (n < 4 || n >5)
      error("ERROR: setPWM() does not have width of 4!");
  }
  switch (type)
  {
    case PCM:
      PCM2PFM();
      PFM2PSSM();
      break;
    case PFM:
      PFM2PSSM();
      break;
    case PSSM:
      break;
    case BEM:
      break;
    default:
      error("setPWM() unrecognized pwm type");
      break;
  }
  setNscore();
  calc_max_score();
}

void PWM::setGC(double gc)
{
  if (input_type == -1)
    this->gc = gc;
  else
  {
    vector<vector<double> > replacement = getPWM(input_type);
    this->gc = gc;
    setPWM(replacement, input_type);
  }
}

void PWM::setPseudo(double pseudo)
{
  if (input_type == -1)
    this->pseudo = pseudo;
  else
  {
    vector<vector<double> > replacement = getPWM(input_type);
    this->pseudo = pseudo;
    setPWM(replacement, input_type);
  }
}

  
vector<vector<double> >& PWM::getPWM()
{
  return mat->getValue();
}

vector<vector<double> > PWM::getPWM(int type)
{
  // work with a copy of the matrix to start
  vector<vector<double> > out = mat->getValue();
  
  switch (type)
  {
    case PCM:
      PSSM2PFM(out);
      PFM2PCM(out);
      break;
    case PFM:
      PSSM2PFM(out);
      break;
    case PSSM:
      calc_max_score();
      break;
    case BEM:
      calc_max_score();
      break;
    default:
      error("getPWM() unrecognized pwm type " + to_string_(type));
      break;
  }
  return out;
}

void PWM::getParameters(param_ptr_vector& p)
{
  if (is_pwm)
  {
    if (mat->isAnnealed())
      p.push_back(mat);
  }
  else if (is_periodic)
  {
    if (period->isAnnealed())
      p.push_back(period);
    if (beta->isAnnealed())
      p.push_back(beta);
  }
  else
    error("undefined type of binding preference");
}

void PWM::getAllParameters(param_ptr_vector& p)
{
  if (is_pwm)
    p.push_back(mat);
  else if (is_periodic)
  {
    p.push_back(period);
    p.push_back(beta);
  }
  else
    error("undefined type of binding preference");
}
  
void PWM::PCM2PFM()
{
  // we have counts. We need to add pseudo count and divide rows
  vector<vector<double> >& matrix = mat->getValue();
  int pwmlen = matrix.size();
  position_counts.resize(pwmlen);
  double gc_adjusted_count;
  
  for (int i=0; i<pwmlen; i++)
  {
    double rowsum = 0;
    for (int j=0; j<4; j++)
    {
      if (j==0 || j==3) // we have an A or T
        gc_adjusted_count = pseudo*(1-gc)/2;
      else             // we have a G or C
        gc_adjusted_count = pseudo*(gc)/2;

      rowsum       += matrix[i][j];
      matrix[i][j] += gc_adjusted_count;
    }
    position_counts[i] = rowsum;
    rowsum += pseudo;

    for (int j=0; j<4; j++)
      matrix[i][j] /= rowsum;
  }
}

void PWM::PFM2PCM(vector<vector<double> >& t)
{
  unsigned int pwmlen = t.size();
  
  if (position_counts.size() != pwmlen)
  {
    warning("Position counts not set. Using default of 100.");
    position_counts.resize(pwmlen);
    for (unsigned int i=0; i<pwmlen; i++)
      position_counts[i] = 100;
  }
  
  double gc_adjusted_count;
  for (unsigned int i=0; i<pwmlen; i++)
  {
    double rowsum = position_counts[i] + pseudo;
    for (int j=0; j<4; j++)
    {
      t[i][j] *= rowsum;
      
      if (j==0 || j==3) // we have an A or T
        gc_adjusted_count = pseudo*(1-gc)/2;
      else             // we have a G or C
        gc_adjusted_count = pseudo*(gc)/2;

      t[i][j] -= gc_adjusted_count;
    }
  }
}

void PWM::PFM2PSSM()
{
  vector<vector<double> >& matrix = mat->getValue();
  int pwmlen = matrix.size();
  double bkgd;
  
  maxscore = 0;
  
  for (int i=0; i<pwmlen; i++)
  {
    double position_max = 0;
    for (int j=0; j<4; j++)
    {
      if (j==0 || j==3)   // we have an A or T
        bkgd = (1-gc)/2;
      else                // we have a G or C
        bkgd = (gc)/2;
      
      matrix[i][j] = log( matrix[i][j]/bkgd );
      position_max = max(position_max, matrix[i][j]);
    }
    maxscore += position_max;
  }
}

void PWM::calc_max_score()
{
  if (is_pwm)
  {
    vector<vector<double> >& matrix = mat->getValue();
    int pwmlen = matrix.size();
    consensus.resize(pwmlen);
    maxscore = 0;
    for (int i=0; i<pwmlen; i++)
    {
      int    pconsensus   = 0;
      double position_max = 0;
      
      for (int j=0; j<4; j++)
      {
        if (matrix[i][j] > position_max)
        {
          position_max = matrix[i][j];
          pconsensus = j;
        }
      }
      consensus[i] = pconsensus;
      maxscore += position_max;
    }
  }
  else if (is_periodic)
  {
    //double pi = 3.14159265358979323846;
    maxscore = 0;
    //ndist = plength/2; // returns floor for middle position
    //mdist = plength-ndist;
    //for (int i = -ndist; i<=mdist; i++)
    //{
    //  double m = 0.25;
    //  double x;
    //  x = 0.25 + beta->getValue() * cos(2*pi*(i/period->getValue()));
    //  m = max(m, x);
    //  m = max(m, (1-x)/2);
    //  x = 0.25 + beta->getValue() * cos(2*pi*(i/period->getValue() + 0.5));
    //  m = max(m, x);
    //  m = max(m, (1-x)/3);
    //  maxscore += log(m/0.25);
    //}
  }
}

/* this will check to make sure frequencies add to 1 */
void PWM::PSSM2PFM(vector<vector<double> >& t)
{
  int pwmlen = t.size();
  double bkgd;
  
  for (int i=0; i<pwmlen; i++)
  {
    double sum = 0;    
    for (int j=0; j<4; j++)
    {
      if (j==0 || j==3)   // we have an A or T
        bkgd = (1-gc)/2;
      else                // we have a G or C
        bkgd = (gc)/2;
      
      t[i][j] = bkgd * exp(t[i][j]);
      sum += t[i][j];
    }
    for (int j=0; j<4; j++)
      t[i][j] /= sum;
  }
}

void PWM::setNscore()
{
  vector<vector<double> >& matrix = mat->getValue();
  int pwmlen = matrix.size();
  for (int i=0; i<pwmlen; i++)
  {
    double minscore = matrix[i][0];
    for (int j=1; j<4; j++)
      minscore = min(minscore, matrix[i][j]);
    matrix[i][4] = minscore;
  }
}

/* returns the threshold that will give this p value. It works... but I have 
no idea how. It was taken from MOODs in BioPerl */
double PWM::pval2score(double p)
{
  vector<vector<double> >& matrix = mat->getValue();
  int n = matrix.size();

  double dp_multi = 100.0;
  
  vector<double> bg;
  bg.push_back( (1-gc)/2 );
  bg.push_back( (  gc)/2 );
  bg.push_back( (  gc)/2 );
  bg.push_back( (1-gc)/2 );

  vector<vector<int> > tmat(4,vector<int>(n,0));

  int maxT = 0;
  int minV = INT_MAX;

  for (int i = 0; i < n; ++i)
  {
      for (int j = 0; j < 4; ++j)
      {
          if (matrix[i][j] > 0.0){
              tmat[j][i] = (int) ( dp_multi * matrix[i][j] + 0.5 );
          }
          else {
              tmat[j][i] = (int) ( dp_multi * matrix[i][j] - 0.5 );
          }
      }
  }

  for (int i = 0; i < n; ++i)
  {
      int max = tmat[0][i];
      int min = max;
      for (int j = 1; j < 4; ++j)
      {
          int v = tmat[j][i];
          if (max < v)
              max = v;
          else if (min > v)
              min = v;
      }
      maxT += max;
      if (minV > min)
          minV = min;
  }

  int R = maxT - n * minV;
  
  vector<double> table0(R + 1, 0.0);
  vector<double> table1(R + 1, 0.0);

  for (int j = 0; j < 4; ++j)
      table0[tmat[j][0] - minV] += bg[j];

  for (int i = 1; i < n; ++i)
  {
      for (int j = 0; j < 4; ++j)
      {
          int s = tmat[j][i] - minV;
          for (int r = s; r <= R; ++r)
              table1[r] += bg[j] * table0[r - s];
      }
      for (int r = 0; r <= R; ++r)
      {
          table0[r] = table1[r];
          table1[r] = 0.0;
      }
  }

  double sum = 0.0;
  
  for (int r = R; r >= 0; --r)
  {
      sum += table0[r];
      if (sum > p)
      {
          return (double) ((r + n * minV + 1) / dp_multi);
      }
  }

	return (double) ((n * minV) / dp_multi);
}

/* a rather inefficient way of converting a score to a pvalue */
double PWM::score2pval(double s)
{
  if (s > maxscore)
    error("score greater than maxscore!");
  
  //calculates log10 p vals up to 10
  int precision = 1000; // number of steps to take
  if (s2p.size() == 0)
  {
    s2p.resize(precision);
    for (int i=1; i<precision; i++)
    {
      double mlog10p = (double) i/(precision/10.0); 
      double pval = pow(10.0,-mlog10p);
      s2p[i] = pval2score(pval);
    }
  }
      
  for (int i=(precision-1); i>0; i--)
  {
    if (s2p[i] < s)
    {
      double mlog10p = (double) i/(precision/10.0); 
      double pval = pow(10.0,-mlog10p);
      return(pval);
    }
  }
  return 1/std::pow(4.0, (int) mat->getValue().size()); // the max pvalue any matrix can have
}   

void PWM::subscore(const vector<int> & s, double * out)
{
  if (is_pwm)
  {
    int i,j,k;
    out[0]=0.0;
    out[1]=0.0;
    vector<vector<double> >& matrix = mat->getValue();
    int pwmlen = matrix.size();
    for (i=0, j=(pwmlen-1); i<pwmlen; i++, j--)
    {
      vector<double>& row = matrix[i]; 
      // forward sequence, i iterates forward
      k = s[i];
      out[0] += row[k];
      
      // reverse sequence, j iterates back over subseq, 3-n give complement
      if (s[j] !=4) 
        k = 3-s[j];
      else 
        k = 4;
      out[1] += row[k];
    }
    out[2] = max(out[0],out[1]);
  }
  else if (is_periodic)
  {
    int i, j;
    int size     = s.size();
    double halfsize = (double) (size - 1) / 2.0; // the position of the dyad, the center of the sequence
      
    out[0] = 0.0;
    out[1] = 0.0;
      
    for (i=1, j=(size-2); i<size; i++, j--)
    {
      // get position with respect to the diad
      double diad_position_f = i - halfsize;
      double diad_position_r = j - halfsize;
      
      // what was the last position
      int first_f  = s[i-1];
      int second_f = s[i];
      
      int first_r  = s[j];
      int second_r = s[j+1];
      
      if (first_r != 4)
        first_r = 3 - first_r;
      if (second_r != 4)
        second_r = 3 - second_r;
      
      out[0] += log(score_dyad(first_f, second_f, diad_position_f)/0.25);
      out[1] += log(score_dyad(first_r, second_r, diad_position_r)/0.25);
    }
    out[2] = max(out[0],out[1]);
  }
  else
    error("no score function for type");
      
}

double PWM::score_dyad(int first, int second, double position)
{
  double x;
  double pi = 3.14159265358979323846;
  if (first == 0)
  {
    x = 0.25 + beta->getValue() * cos(2 * pi * (position/period->getValue()));
    if (second == 0)
      return x;
    else
      return (1-x)/3;
  }
  else if (first == 3)
  {
    x = 0.25 + beta->getValue() * cos(2 * pi * (position/period->getValue()));
    if (second == 0 || second == 3)
      return x;
    else
      return (1-2*x)/2;
  }
  else if (first == 2)
  {
    x = 0.25 + beta->getValue() * cos(2*pi*(position/period->getValue() + 0.5));
    if (second == 1)
      return x;
    else
      return (1-x)/3;
  }
  else
    return 0.25;
}
      
  

void PWM::score(const vector<int>& s, TFscore &t)
{
  /* I do something a little unorthodox here. I want to control for boundary
  effects, so I pad the beginning and end of the sequence with Ns and score
  using the boundary behavior. Unlike the transc code, I report the score for the
  middle of the binding site rather than the m position. This means my output is
  of the same length as the sequence and it is the same length for every factor */
  
  int i, j, k;
  int mdist, ndist;
  int start, end;
  int len;
  int pwmlen;
  
  if (is_pwm)
  {
    vector<vector<double> >& matrix = mat->getValue();
    pwmlen = matrix.size();
  }
  else if (is_periodic)
    pwmlen = plength;
  else
    error("unrecognized pwm type in score");
    
    
  vector<int> sub(pwmlen);
  ndist = pwmlen/2; // returns floor for middle position
  mdist = pwmlen-ndist;
  
  double * out = new double[3];
  len = s.size();

  t.fscore.resize(len);
  t.rscore.resize(len);
  t.mscore.resize(len);
  //t.tfname = tfname;
  
  for (i=0; i<len; i++)
  {
    start = i - mdist;
    end   = i + ndist;
    // j = index in s
    // j = index in pwm
    for (j=start, k=0; j<end; j++, k++)
    {
      if (j<0 || j>=len)
      {
        sub[k]=4;
      } else 
      {
        sub[k]=s[j];
      }
    }
    subscore(sub,out);
    t.fscore[i] = out[0];
    t.rscore[i] = out[1];
    t.mscore[i] = out[2];
  } 
  delete[] out;
}

//size_t PWM::getSize()
//{
//  int n = mat.size()*4 + 1;
//  return n*sizeof(int);
//}
//
//void PWM::serialize(void *buf) const
//{
//  int index = 0;
//  int n = mat.size();
//  int* dest = static_cast<int *>(buf);
//  dest[index++] = n;
//  for (int i=0; i<n; i++)
//  {
//    for (int j=0; j<4; j++)
//      dest[index++] = mat[i][j];
//  }
//}
//      
//void PWM::deserialize(void const *buf)
//{
//  int index = 0;
//  int const * from = static_cast<int const *>(buf);
//  int n = from[index++];
//  for (int i=0; i<n; i++)
//  {
//    for (int j=0; j<4; j++)
//      mat[i][j] = from[index++];
//  }
//  setNscore();
//  calc_max_score();
//}


//void PWM::print(ostream& os, int precision)
//{
//  int p = precision;
//  int w = p + 7;
//  vector<vector<do
//  int pwmlen = mat.size();
//  vector<char> bases;
//  bases.push_back('A');
//  bases.push_back('C');
//  bases.push_back('G');
//  bases.push_back('T');
//
//  for (int i=0; i<4; i++)
//  {
//    os << setw(w) << setprecision(p) << bases[i];
//    for (int j=0; j<pwmlen; j++)
//      os << setw(w) << setprecision(p) << mat[j][i];
//    os << endl;
//  }
//}
