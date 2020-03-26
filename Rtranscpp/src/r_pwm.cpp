#include "r_pwm.h"

using namespace Rcpp;

PWMPtr::PWMPtr(NumericMatrix mat, string type)
: pwm(new PWM)
{
  set_pwm(mat, type);
}

PWMPtr::PWMPtr(NumericMatrix mat, string type, double gc, double pseudo)
: pwm(new PWM)
{
  set_pwm(mat, type);
  pwm->setGC(gc);
  pwm->setPseudo(pseudo);
}
  
  
NumericMatrix PWMPtr::get_pwm(string type)
{
  PWM& p = *pwm;
  vector<vector<double> > mat;
  
  cerr << "getting pwm" << endl;
  if (type == string("PCM"))
    mat = p.getPWM(PCM);
  else if (type == string("PFM"))
    mat = p.getPWM(PFM);
  else if (type == string("PSSM"))
    mat = p.getPWM(PSSM);

  int pwmlen = mat.size();
  
  NumericMatrix m(4, pwmlen);
  for (int i=0; i<pwmlen; i++)
  {
    for (int j=0; j<4; j++)
      m(j,i) = mat[i][j];
  }
  return m;
}

void PWMPtr::set_pwm(NumericMatrix input, string type)
{
  PWM& p = *pwm;
  vector<vector<double> > mat;
  
  int nrow = input.nrow();
  int ncol = input.ncol();
  
  if (nrow != 4)
    error("Input matrix must have 4 rows (A C G T)");
  
  mat.resize(ncol);
  for (int i=0; i<ncol; i++)
  {
    for (int j=0; j<4; j++)
      mat[i].push_back(input(j,i));
  }
  
  if (type == string("PCM"))
    p.setPWM(mat, PCM);
  else if (type == string("PFM"))
    p.setPWM(mat, PFM);
  else if (type == string("PSSM"))
    p.setPWM(mat, PSSM);
}     

List PWMPtr::score(CharacterVector seq)
{
  vector<char> char_seq;
  
  /* the conversion to character vector really doesnt seem to care if the
  input is a string or character vector so this accepts either input */
  for (int i=0; i<seq.size(); i++)
    char_seq.push_back(Rcpp::as<char>(seq[i]));
  
  TFscore output;
  vector<int> int_seq = char2int(char_seq);
  pwm->score(int_seq, output);
  
  List out(4);
  out[0] = output.fscore;
  out[1] = output.rscore;
  out[2] = output.mscore;
  out[3] = pwm->getMaxScore();
  
  vector<string> names;
  names.push_back("f");
  names.push_back("r");
  names.push_back("m");
  names.push_back("maxscore");
  out.attr("names") = names;
  
  return out;
}



RCPP_MODULE(mod_pwm)
{
  class_<PWMPtr>("PWM")
  
  .constructor()
  .constructor<NumericMatrix,string>()
  .constructor<NumericMatrix,string,double,double>()
  
  .property("source", &PWMPtr::get_source , &PWMPtr::set_source )
  .property("gc",     &PWMPtr::get_gc     , &PWMPtr::set_gc     ) 
  .property("pseudo", &PWMPtr::get_pseudo , &PWMPtr::set_pseudo )

  .property("length",    &PWMPtr::get_length)
  .property("max_score", &PWMPtr::get_max)
  .property("pcm",       &PWMPtr::get_pcm)
  .property("pfm",       &PWMPtr::get_pfm)
  .property("pssm",      &PWMPtr::get_pssm)
  
  .method("get_pwm", &PWMPtr::get_pwm)
  .method("set_pwm", &PWMPtr::set_pwm)
  
  .method("pval2score", &PWMPtr::pval2score)
  .method("score2pval", &PWMPtr::score2pval)
  
  .method("score", &PWMPtr::score)
  
  ;
}

