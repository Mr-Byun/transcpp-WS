/*********************************************************************************
*                                                                                *
*     bindingsite.cpp                                                            *
*                                                                                *
*     Contains the structure definition for a binding site.                      *
*     Pretty dull for now                                                        *
*                                                                                *
*********************************************************************************/

#include "bindingsite.h"
#include <iostream>
#include <iomanip>
#include "TF.h"

void printSiteHeader(ostream& os)
{
  int w = 10;
  os << setprecision(3)
     << setw(w) << "name"
     << setw(w) << "start"
     << setw(w) << "end"
     << setw(w) << "score"
     << setw(w) << "K"
     << setw(w) << "K*kmax"
     << setw(w) << "strand";
}

void printSite(const BindingSite& b, ostream& os)
{
  int w = 10;
  os << setprecision(3)
     << setw(w) << b.tf->getName()
     << setw(w) << b.m 
     << setw(w) << b.n 
     << setw(w) << b.score 
     << setw(w) << b.K_exp_part
     << setw(w) << b.K_exp_part_times_kmax;
  if (b.orientation == 'F')
    os << setw(w) << "+";
  else 
    os << setw(w) << "-";
}

void printSite(BindingSite* b, ostream& os) { printSite(*b, os); }

/*
bool compareBindingSiteRight(BindingSite* a, BindingSite* b)
{
  return(a->n < b->n);
}
*/

bool compareBindingSiteRight(const BindingSite* a, const BindingSite* b)
{
  int A = a->n;
  int B = b->n;
    
  if (A < B)
    return true;
  else if (A > B)
    return false;
  else if (a==b)
    return false;
  else
  {
    TF& tfa = *(a->tf);
    TF& tfb = *(b->tf);
    int idxa = tfa.getIndex();
    int idxb = tfb.getIndex();
    
    if (idxa < idxb)
      return true;
    else if (idxa > idxb)
      return false;
    else if (a->orientation == 'F' && b->orientation == 'R')
      return true;
    else if (a->orientation == 'R' && b->orientation == 'F')
      return false;
    else
    {
      printSiteHeader(cerr);
      cerr << endl;
      printSite(*a, cerr);
      cerr << endl;
      printSite(*b, cerr);
      cerr << endl;
      error("sites should never be equal!");
      return false;
    }
  }
}

/*
bool compareBindingSiteLeft(BindingSite* a, BindingSite* b)
{
  return(a->m > b->m);
}
*/

bool compareBindingSiteLeft(BindingSite* a, BindingSite* b)
{
  int A = a->m;
  int B = b->m;
  
  if (A > B)
    return true;
  else if (A < B)
    return false;
  else if (a==b)
    return false;
  else
  {
    TF& tfa = *(a->tf);
    TF& tfb = *(b->tf);
    int idxa = tfa.getIndex();
    int idxb = tfb.getIndex();
    if (idxa < idxb)
      return true;
    else if (idxa > idxb)
      return false;
    else if (a->orientation == 'F' && b->orientation == 'R')
      return true;
    else if (a->orientation == 'R' && b->orientation == 'F')
      return false;
    else
    {
      error("sites should never be equal!");
      return false;
    }
  }
}

bool overlaps(BindingSite& site1, BindingSite& site2)
{
  return ( site1.m < site2.n && site2.m < site1.n);
}

bool bad_overlap_function(BindingSite& site1, BindingSite& site2)
{
  return ( site1.m <= site2.n && site2.m <= site1.n);
}
