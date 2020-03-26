#include <Rcpp.h>
#include <boost/shared_ptr.hpp>
#include "organism.h"

#ifndef R_DATATABLE_H
#define R_DATATABLE_H

using namespace Rcpp;

class DataTablePtr
{
private:
  boost::shared_ptr<DataTable<double> > table;
  
public:
  DataTablePtr();
  DataTablePtr(string fname, string tname);
  DataTablePtr(string fname, string tname, string section);
  DataTablePtr(table_ptr table) { this->table = table; } // cannot be exposed!
  
  List data_frame();
  void set(DataFrame x, string row_type, string col_type);
};






















#endif
