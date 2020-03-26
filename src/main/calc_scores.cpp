/*********************************************************************************
*                                                                                *
*     calc_scores.cpp                                                            *
*                                                                                *
*     For cases where we dont know concentrations, calculate binding scores      *
*     on DNA and skip all other calculations                                     *
*                                                                                *
*********************************************************************************/

#include "pwm.h"
#include "fasta.h"
#include "gene.h"
#include "TF.h"
#include "parameter.h"
#include "utils.h"

#include <unistd.h>
#include <getopt.h>

static const char *optString = "f:p:";

static const struct option longOpts[] = {
    { "help",        no_argument,       NULL, 'h' },
    { "fasta",       required_argument, NULL, 'f' },
    { "pwms",        required_argument, NULL, 'p' },
    { "p-value",     no_argument,       NULL,  0  },
    { "gene",        required_argument, NULL,  0  },
    { "tf",          required_argument, NULL,  0  },
    { "tfnames",     no_argument,       NULL,  0  },
    { "gnames",      no_argument,       NULL,  0  },
    { "maxscores",   no_argument,       NULL,  0  },
    { "forward",     no_argument,       NULL,  0  },
    { "reverse",     no_argument,       NULL,  0  },
    { "normalize",   no_argument,       NULL,  0  },
    { 0, 0, 0, 0}
};

void display_usage()
{
  cerr << endl << "\t Usage" << endl << endl
       << "\t calc_scores -f [fasta-file] -p [pwm-file] [options]" << endl << endl
       << "\t Options" << endl
       << "\t --help    [-h]   print this message" << endl
       << "\t --fasta   [-f]   a fasta formatted file containing input sequences" << endl
       << "\t --pwms    [-p]   the file containing pwms to use" << endl
       << "\t --p-value        calculate -log10 p-values (slow)" << endl
       << "\t --gene           calculate for specified gene" << endl
       << "\t --tf             calculate for specified pwm" << endl
       << "\t --tfnames        print tf names and exit" << endl
       << "\t --gnames         print gene names and exit" << endl
       << "\t --forward        prints for forward strand" << endl
       << "\t --reverse        prints for reverse strand" << endl
       << "\t --maxscores      prints the max score for each tf" << endl
       << "\t --normalize      normalize scores to 0=min, 1=max" << endl << endl;

  exit(1);
}

int main(int argc, char* argv[])
{
  int opt = 0;
  int longIndex = 0;
  
  bool   pval          = false;
  bool   normalize     = false;
  bool   print_tfnames = false;
  bool   print_gnames  = false;
  bool   forward       = false;
  bool   reverse       = false;
  bool   maxscores     = false;
  
  string fasta_file("");
  string pwm_file("");
  string gname("");
  string tfname("");
  
  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while(opt != -1)
  {
    switch (opt)
    {
      case 'h':
        display_usage();
        break;
      case 'f':
        fasta_file = optarg;
        break;
      case 'p':
        pwm_file   = optarg;
        break;
      case 0:
        if (longOpts[longIndex].name == "p-value")
          pval      = true;
        else if (longOpts[longIndex].name == "normalize")
          normalize = true;
        else if (longOpts[longIndex].name == "gene")
          gname = optarg;
        else if (longOpts[longIndex].name == "tf")
          tfname = optarg;
        else if (longOpts[longIndex].name == "tfnames")
          print_tfnames = true;
        else if (longOpts[longIndex].name == "gnames")
          print_gnames = true;
        else if (longOpts[longIndex].name == "forward")
          forward = true;
        else if (longOpts[longIndex].name == "reverse")
          reverse = true;
        else if (longOpts[longIndex].name == "maxscores")
          maxscores = true;

        break;
      case '?':
        display_usage();
        break;
      default:
        display_usage();
        break;
    }
    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  }
  
  if (fasta_file == string("") || 
      pwm_file   == string(""))
  {
    cerr << "No fasta or pwm file specified" << endl;
    display_usage();
  }
  
  // read the files
  PwmContainer pwms(pwm_file);
  Fasta fasta(fasta_file);
  
  // create genes from fasta file
  GeneContainer genes;
  vector<string>& gnames = fasta.getNames();
  int ngenes = gnames.size();
  
  if (print_gnames)
  {
    cout << "$gene_names" << endl;
    for (int i=0; i<ngenes; i++)
      cout << gnames[i] << endl;
    exit(1);
  }
  
  for (int i=0; i<ngenes; i++)
  {
    gene_ptr gene(new Gene());
    gene->setName(gnames[i]);
    gene->setSequence(fasta.getSeq(gnames[i]));
    genes.add(gene);
  }
  
  // create TFs from pwm file
  TFContainer tfs;
  vector<string>& tfnames = pwms.getNames();
  int ntfs = tfnames.size();
  
  if (print_tfnames)
  {
    cout << "$tf_names" << endl;
    for (int i=0; i<ntfs; i++)
      cout << tfnames[i] << endl;
    exit(1);
  }
  
  for (int i=0; i<ntfs; i++)
  {
    tf_ptr tf(new TF());
    tf->setName(tfnames[i]);
    tf->setPwm(pwms.getPwm(i).getLogOddPwm(1,0.407),0.407);
    tfs.add(tf);
  }
  
  if (maxscores)
  {
    int namesize = 0;
    for (int i=0; i<ntfs; i++)
      namesize = max(namesize, (int) tfnames[i].size());
  
    int precision = 4;
    int w1 = namesize+1;
    int w2 = precision + 5;
  
    cout << "$maxscores" << endl;
    for (int i=0; i<ntfs; i++)
    {
      TF& tf = tfs.getTF(i);
      cout << setw(w1) << tf.getName();
      cout << setw(w2) << setprecision(precision) << tf.getMaxScore() << endl;
    }
    exit(1);
  }
      
  if (gname==string(""))
  {
    cerr << "no gene specified" << endl;
    exit(1);
  }
  
  if (tfname==string(""))
  {
    cerr << "no tf specified" << endl;
    exit(1);
  }
  
  if (forward && reverse)
  {
    cerr << "cannot print both forward and reverse strand!" << endl;
    exit(1);
  }
  
  Gene& gene = genes.getGene(gname);
  TF&   tf   = tfs.getTF(tfname);
  
  TFscore tfscore;
  tf.score(gene.getSequence(), tfscore);
  int length = tfscore.mscore.size();
  
  vector<double> scores;
  
  if (forward)
    scores = tfscore.fscore;
  else if (reverse)
    scores = tfscore.rscore;
  else
    scores = tfscore.mscore;
  
   cout << "$scores" << endl;
   if (pval)
   {
     for (int i=0; i<length; i++)
       scores[i] = tf.score2pval(scores[i]);
   }
   for (int i=0; i<length; i++)
     cout << scores[i] << endl;

  return 0;
    
    
  /*
  // do the clalculations
  for (int i=0; i<ngenes; i++)
  {
    Gene& gene = genes.getGene(i);
    map<TF*, TFscore>& gscores = scores[&gene];
    //map<TF*, vector<double > >& gpvals  = pvals[&gene];
    for (int j=0; j<ntfs; j++)
    {
      TF& tf = tfs.getTF(j);
      TFscore& tfscore = gscores[&tf];
      //vector<double>& tfpvals = gpvals[&tf];
      tf.score(gene.getSequence(), tfscore);
      int length = tfscore.mscore.size();
      //tfpvals.resize(length);
      //for (int k=0; k<length; k++)
      //  tfpvals[k] = tf.score2pval(tfscore.mscore[k]);
    }
  }
  
  // find the width of names
  int namesize = 0;
  for (int i=0; i<ntfs; i++)
    namesize = max(namesize, (int) tfnames[i].size());
  
  int precision = 3;
  int w = max(namesize+1, precision+5);
    
    
  for (int i=0; i<1; i++)
  {
    Gene& gene = genes.getGene(i);
    int length = gene.length();
    // print the header
    cout << '$' << gene.getName() << '_' << "scores" << endl;
    cout << setw(w) << "bp";
    for (int j=0; j<ntfs; j++)
      cout << setw(w) << tfnames[j];
    cout << endl;
    // print the data
    for (int j=0; j<length; j++)
    {
      cout << setw(w) << j;
      for (int k=0; k<ntfs; k++)
      {
        TF& tf = tfs.getTF(k);
        cout << setw(w) << scores[&gene][&tf].mscore[j];
      }
      cout << endl;
    }
  }*/
}
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
          
        
