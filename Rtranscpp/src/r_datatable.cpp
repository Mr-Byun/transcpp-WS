#include <Rcpp.h>
#include "r_datatable.h"

using namespace Rcpp;

DataTablePtr::DataTablePtr() :
  table(boost::shared_ptr<DataTable<double> >(new DataTable<double>())) 
{}


DataTablePtr::DataTablePtr(string fname, string tname) :
  table(boost::shared_ptr<DataTable<double> >(new DataTable<double>()))
{
  ptree pt;
  read_xml(fname, pt, boost::property_tree::xml_parser::trim_whitespace);
  ptree& root_node  = pt.get_child("Root");
  ptree& input_node = root_node.get_child("Output");

  table->read(input_node, tname);
}

DataTablePtr::DataTablePtr(string fname, string tname, string section) :
  table(boost::shared_ptr<DataTable<double> >(new DataTable<double>()))
{
  ptree pt;
  read_xml(fname, pt, boost::property_tree::xml_parser::trim_whitespace);
  ptree& root_node  = pt.get_child("Root");
  ptree& input_node = root_node.get_child(section);

  table->read(input_node, tname);
}

List DataTablePtr::data_frame()
{
  vector<string>& row_names = table->getRowNames();
  vector<string>& col_names = table->getColNames();
  
  int nrow = row_names.size();
  int ncol = col_names.size();
  
  vector<double> t_data;
  t_data.resize(nrow);
  
  List out(ncol);
  
  for (int i=0; i<ncol; i++)
  {
    string& col_name = col_names[i];
    vector<double*>& col = table->getCol(col_name);
    for (int j=0; j<nrow; j++)
      t_data[j] = *(col[j]);
    NumericVector new_col = wrap(t_data);
    out[i] = new_col;
  }
  out.attr("class")     = string("data.frame");
  out.attr("row.names") = row_names;
  out.attr("names")     = col_names;
  return out;
}

void DataTablePtr::set(DataFrame x, string row_type, string col_type)
{
  vector<string> row_names = x.attr("row.names");
  vector<string> col_names = x.attr("names");
  
  int nrow = row_names.size();
  int ncol = col_names.size();
  
  vector< vector<double> > data(nrow, vector<double>(ncol));
  for (int i=0; i<ncol; i++)
  {
    NumericVector col = x[i];
    for (int j=0; j<nrow; j++)
    {
      data[j][i] = col[j];
    }
  }
  table->set(data, row_names, col_names, row_type, col_type);
    
}

  
RCPP_MODULE(mod_datatable)
{
  class_<DataTablePtr>("DataTable")
  
  .constructor()
  .constructor<string, string>()
  .constructor<string, string, string>()
  
  .property("data_frame", &DataTablePtr::data_frame)
  
  .method("set", &DataTablePtr::set)
  ;
}
  
  
