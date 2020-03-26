#include "mex.h"
#include "../src/pwm.h"
#include "class_handle.hpp"
#include "../src/organism.h"

void create(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Check parameters
  if (nlhs != 1)
    mexErrMsgTxt("New: One output expected.");
  if (nrhs < 2)
    mexErrMsgTxt("New: One input expected.");
  // get the string filename as input and convert to ptree

  string fname(mxArrayToString(prhs[1]));

  ptree pt;
  read_xml(fname, pt, boost::property_tree::xml_parser::trim_whitespace);

  ptree& root_node   = pt.get_child("Root");
  ptree& mode_node   = root_node.get_child("Mode");

  mode_ptr mode(new Mode(fname, mode_node));

  // check if there is a second argument
  string section_name;
  if (nrhs > 2)
    section_name = mxArrayToString(prhs[2]);
  else
    section_name = "Output";

  ptree& section_node = root_node.get_child(section_name);

  // Return a handle to a new C++ instance
  mexPrintf("Reading section %s from file %s \n", section_name.c_str(), fname.c_str());

  plhs[0] = convertPtr2Mat<Organism>(new Organism(section_node, mode));
  return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Get the command string
  char cmd[64];
  if (nrhs < 1 || mxGetString(prhs[0], cmd, sizeof(cmd)))
    mexErrMsgTxt("First input should be a command string less than 64 characters long.");

  // New
  if (!strcmp("new", cmd))
  {
    create(nlhs, plhs, nrhs, prhs);
    return;
  }

  // Check there is a second input, which should be the class instance handle
  if (nrhs < 2)
    mexErrMsgTxt("Second input should be a class instance handle.");

  // Delete
  if (!strcmp("delete", cmd))
  {
    // Destroy the C++ object
    destroyObject<Organism>(prhs[1]);
    // Warn if other commands were ignored
    if (nlhs != 0 || nrhs != 2)
      mexWarnMsgTxt("Delete: Unexpected arguments ignored.");
    return;
  }

  // Get the class instance pointer from the second input
  Organism *Organism_instance = convertMat2Ptr<Organism>(prhs[1]);

  // Call the various class methods

  // Test
  if (!strcmp("test", cmd))
  {
    // Check parameters
    if (nlhs < 0 || nrhs < 2)
      mexErrMsgTxt("Test: Unexpected arguments.");
    // Call the method
    // Organism_instance->test();
    plhs[0] = mxCreateDoubleScalar(Organism_instance->test());

    return;
  }

  if (!strcmp("getNNuc", cmd))
  {
    // Check parameters
    if (nlhs < 0 || nrhs < 2)
      mexErrMsgTxt("getNNuc: Unexpected arguments.");
    // Call the method
    // Organism_instance->test();
    plhs[0] = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
    int* output = (int*) mxGetData(plhs[0]);
    output[0] = Organism_instance->getNNuc();

    return;
  }

  if (!strcmp("getNGenes", cmd))
  {
    // Check parameters
    if (nlhs < 0 || nrhs < 2)
      mexErrMsgTxt("getNGenes: Unexpected arguments.");
    // Call the method
    // Organism_instance->test();
    plhs[0] = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
    int* output = (int*) mxGetData(plhs[0]);
    output[0] = Organism_instance->getNGenes();

    return;
  }

  if (!strcmp("getNParams", cmd))
  {
    // Check parameters
    if (nlhs < 0 || nrhs < 2)
      mexErrMsgTxt("getNParams: Unexpected arguments.");
    // Call the method
    // Organism_instance->test();
    plhs[0] = mxCreateNumericMatrix(1,1,mxUINT32_CLASS,mxREAL);
    int* output = (int*) mxGetData(plhs[0]);
    output[0] = Organism_instance->getDimension();

    return;
  }

  if (!strcmp("getTotalScore", cmd))
  {
    // Check parameters
    if (nlhs < 0 || nrhs < 2)
      mexErrMsgTxt("getTotalScore: Unexpected arguments.");
    // Call the method
    // Organism_instance->test();
    plhs[0] = mxCreateDoubleScalar(Organism_instance->getTotalScore());
    return;
  }

  if (!strcmp("getParameterName", cmd))
  {
    // Check parameters
    if (nlhs < 0 || nrhs < 3)
      mexErrMsgTxt("getParameterName: Unexpected arguments.");
    // Call the method
    double* dinput = (double*) mxGetData(prhs[2]);
    int     idx    = (int) dinput[0];
    idx = idx - 1;

    int nparams = Organism_instance->getDimension();
    if (idx >= nparams || idx < 0)
    {
      stringstream err;
      err << "getParameterName: Parameter index " << idx << " out of bounds.";
      mexErrMsgTxt(err.str().c_str());
    }
    string pname   = Organism_instance->getParamName(idx);
    plhs[0] = mxCreateString(pname.c_str());
    return;
  }

  if (!strcmp("getParameterValue", cmd))
  {
    // Check parameters
    if (nlhs < 0 || nrhs < 3)
      mexErrMsgTxt("getParameterValue: Unexpected arguments.");
    // Call the method
    double* dinput = (double*) mxGetData(prhs[2]);
    int     idx    = (int) dinput[0];
    idx--;
    int nparams = Organism_instance->getDimension();
    if (idx >= nparams || idx < 0)
      mexErrMsgTxt("getParameterValue: Parameter index out of bounds.");
    ParameterInterface* param = (Organism_instance->getParam(idx)).get();
    if (param->getType() == string("double"))
    {
      Parameter<double>* p = dynamic_cast<Parameter<double>* >(param);
      plhs[0] = mxCreateDoubleScalar(p->getValue());
    }
    else if (param->getType() == string("Sequence"))
    {
     Parameter<Sequence>* p = dynamic_cast<Parameter<Sequence>* >(param);
      plhs[0] = mxCreateString(p->getValue().getSequenceString().c_str());
    }
    else if (param->getType() == string("PWM"))
    {
      vector<vector<double> > pwm;
      Parameter<PWM>* p = dynamic_cast<Parameter<PWM>* >(param);
      // check number of input params given to see if type was speficied
      if (nrhs > 3)
      {
        string type(mxArrayToString(prhs[3]));
        mexPrintf("%s\n",type.c_str());
        if (type == "PCM")
          pwm = p->getValue().getPWM(0);
        else if (type == "PFM")
          pwm = p->getValue().getPWM(1);
        else if (type == "PSSM")
          pwm = p->getValue().getPWM(2);
        else if (type == "BEM")
          pwm = p->getValue().getPWM(3);
        else
          mexErrMsgTxt("getParameterValue: Unrecognized pwm format requested (use 'PCM', 'PFM', 'PSSM', or 'BEM')");
      }
      else
        pwm = p->getValue().getPWM(3);

      int pwmlen = pwm.size();
      plhs[0] = mxCreateDoubleMatrix((mwSize) 4, (mwSize) pwmlen, mxREAL);
      double* out = mxGetPr(plhs[0]);
      for (int i=0; i<pwmlen; i++)
      {
        vector<double>& row = pwm[i];
        for (int j=0; j<4; j++)
        {
          out[j+i*4] = row[j];
        }
      }
    }
    else
    {
      mexErrMsgTxt("getParameterValue: Unknown parameter type.");
    }
    return;
  }

  if (!strcmp("setParameterValue", cmd))
  {
    // Check parameters
    if (nlhs < 0 || nrhs < 4)
      mexErrMsgTxt("setParameterValue: Unexpected arguments.");
    // Call the method
    double* dinput = (double*) mxGetData(prhs[2]);
    int     idx    = (int) dinput[0];
    idx--;
    int nparams = Organism_instance->getDimension();
    if (idx >= nparams || idx < 0)
      mexErrMsgTxt("setParameterValue: Parameter index out of bounds.");
    ParameterInterface* param = (Organism_instance->getParam(idx)).get();
    if (param->getType() == string("double"))
    {
      Parameter<double>* p = dynamic_cast<Parameter<double>* >(param);
      dinput = (double*) mxGetData(prhs[3]);
      p->set(*dinput);
      Organism_instance->move(idx);
    }
    else if (param->getType() == string("Sequence"))
    {
      Parameter<Sequence>* p = dynamic_cast<Parameter<Sequence>* >(param);
      string input = mxArrayToString(prhs[3]);
      p->set(input);
      Organism_instance->move(idx);
    }
    else if (param->getType() == string("PWM"))
    {
      Parameter<PWM>* p = dynamic_cast<Parameter<PWM>* >(param);
      int pwmlen = mxGetN(prhs[3]);
      int pwmwid = mxGetM(prhs[3]);
      double* in = (double*) mxGetData(prhs[3]);
      if (pwmwid != 4)
        mexErrMsgTxt("setParameterValue: PWM width must be 4");
      vector<vector<double> > pwm;
      for (int i=0; i<pwmlen; i++)
      {
        vector<double> row;
        for (int j=0; j<4; j++)
          row.push_back(in[j+i*4]);
        pwm.push_back(row);
      }
      PWM& pwm_param = p->getValue();
      if (nrhs > 4)
      {
        string type(mxArrayToString(prhs[4]));
        if (type == "PCM")
          pwm_param.setPWM(pwm, 0);
        else if (type == "PFM")
          pwm_param.setPWM(pwm, 1);
        else if (type == "PSSM")
          pwm_param.setPWM(pwm, 2);
        else
          mexErrMsgTxt("getParameterValue: Unrecognized pwm format requested (use 'PCM', 'PFM', or 'PSSM')");
      }
      else
        pwm_param.setPWM(pwm, 2);

      Organism_instance->move(idx);
    }
    else
    {
      mexErrMsgTxt("setParameterValue: Unknown parameter type.");
    }
    return;
  }

  if (!strcmp("getN", cmd))
  {
    int nnuc   = Organism_instance->getNNuc();
    int ngenes = Organism_instance->getNGenes();
    plhs[0] = mxCreateDoubleMatrix((mwSize) nnuc, (mwSize) ngenes, mxREAL);
    double* out = mxGetPr(plhs[0]);
    for (int gene=0; gene<ngenes; gene++)
    {
      vector<double>& Ns = Organism_instance->getN(gene);
      for (int nuc=0; nuc<nnuc; nuc++)
      {
        out[nuc+gene*nnuc] = Ns[nuc];
      }
    }
    return;
  }
  
  if (!strcmp("getR", cmd))
  {
    int nnuc   = Organism_instance->getNNuc();
    int ngenes = Organism_instance->getNGenes();
    plhs[0] = mxCreateDoubleMatrix((mwSize) nnuc, (mwSize) ngenes, mxREAL);
    double* out = mxGetPr(plhs[0]);
    for (int gene=0; gene<ngenes; gene++)
    {
      vector<double>& Rs = Organism_instance->getR(gene);
      for (int nuc=0; nuc<nnuc; nuc++)
      {
        out[nuc+gene*nnuc] = Rs[nuc];
      }
    }
    return;
  }
  
  if (!strcmp("getBindings", cmd))
  {
    
  }

  // Got here, so command not recognized
  mexErrMsgTxt("Command not recognized.");
}



















