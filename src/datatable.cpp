/*********************************************************************************
*                                                                                *
*     datatable.cpp                                                              *
*                                                                                *
*     This contains a simple data-table class. It holds a 2D table and provides  *
*     methods for access and modification. Each row and column has names used    *
*     for access. This also contains methods for reading and writing tables to   *
*     fixed width tables or xml tables. This uses a map implementation with      *
*     string access, so for efficiency inner loop calculations should act on     *
*     pointers to the data contained here rather than calling the get functions  *
*                                                                                *
*********************************************************************************/

#include "datatable.h"
#include "utils.h"

#define foreach_ BOOST_FOREACH

#include <boost/foreach.hpp>
#include <cstdlib>
#include <algorithm>
#include <cstdlib>
#include <boost/version.hpp>

int myrandom (int i) { return std::rand()%i;}

/******************************   DataTable    **********************************/

/*    Constructors    */

template< typename T >
DataTable<T>::DataTable() {}

template< typename T >
DataTable<T>::DataTable(ptree& pt, string node_name)
{
  read(pt, node_name);
}

template< typename T >
DataTable<T>::DataTable(string fname, string node_name, string section)
{
  ptree pt;
  read_xml(fname, pt, boost::property_tree::xml_parser::trim_whitespace);
  ptree& root_node  = pt.get_child("Root");
  ptree& input_node = root_node.get_child(section);
  read(input_node, node_name);
}

/*    Setters   */







/*    Getters   */

template< typename T >
T& DataTable<T>::getDataPoint(string type1, string id1, string type2, string id2) 
{
  if (type1 == row_datatype && type2==col_datatype)
    return data[id1][id2];
  else if (type1 == col_datatype && type2==row_datatype)
    return data[id2][id1];
  else
  {
    stringstream err;
    err << "ERROR: Data table does not contain data types " << type1 << " and " << type2 << endl;
    err << "       Types are " << row_datatype << " and " << col_datatype << endl;
    error(err.str());
    return data[id1][id2]; // you will never get here!
  }
}

template< typename T >
int DataTable<T>::getN(string n)
{
  if (n == row_datatype)
    return row_names.size();
  else if (n == col_datatype)
    return col_names.size();
  else
  {
    stringstream err;
    err << "ERROR: Data table does not contain data type " << n << endl;
    err << "       Types are " << row_datatype << " and " << col_datatype << endl;
    error(err.str());
    return row_names.size(); // you will never get here!
  }
}

template< typename T >
vector<string>& DataTable<T>::getNames(string n)
{
  if (n == row_datatype)
    return row_names;
  else if (n == col_datatype)
    return col_names;
  else
  {
    stringstream err;
    err << "ERROR: Data table does not contain data type " << n << endl;
    err << "       Types are " << row_datatype << " and " << col_datatype << endl;
    error(err.str());
    return row_names; // you will never get here!
  }
}
    
template< typename T >
vector<T*>& DataTable<T>::getRow(string row)
{
  return row_data[row];
}

template< typename T >
vector<T*>& DataTable<T>::getCol(string col)
{
  return col_data[col];
}

  
/*    Methods   */

template< typename T >
bool DataTable<T>::hasRowName(string& n)
{
  int nrows = row_names.size();
  for (int i=0; i<nrows; i++)
  {
    if (row_names[i] == n)
      return true;
  }
  return false;
}

template< typename T >
bool DataTable<T>::hasColName(string& n)
{
  int ncols = col_names.size();
  for (int i=0; i<ncols; i++)
  {
    if (col_names[i] == n)
      return true;
  }
  return false;
}

template< typename T >
void DataTable<T>::fillNaN()
{
  int nrows = row_names.size();
  int ncols = col_names.size();
  for (int i=0; i<nrows; i++)
  {
    string& rowname = row_names[i];
    map<string, T>& row = data[rowname];
    for (int j=0; j<ncols; j++)
    {
      string& colname = col_names[j];
      if (row.find(colname) == row.end())
        row[colname] = numeric_limits<T>::quiet_NaN();
    }
  }
}

template< typename T >
void DataTable<T>::vectorize()
{
  int nrows = row_names.size();
  int ncols = col_names.size();
  
  row_data.clear();
  col_data.clear();
  
  for (int i=0; i<nrows; i++)
  {
    string& rowname = row_names[i];
    map<string, T>& row = data[rowname];
    for (int j=0; j<ncols; j++)
    {
      string& colname = col_names[j];
      row_data[rowname].push_back(&row[colname]);
    }
  }
  
  for (int i=0; i<ncols; i++)
  {
    string& colname = col_names[i];
    for (int j=0; j<nrows; j++)
    {
      string& rowname = row_names[j];
      col_data[colname].push_back(&data[rowname][colname]);
    }
  }
}   

template< typename T >
void DataTable<T>::permute(string& n, int precision)
{
  srand(getpid());
  
  int nrows = row_names.size();
  int ncols = col_names.size();
  
  map<string, map<string, T> > new_data;
  vector<string> new_row_names = row_names;
  vector<string> new_col_names = col_names;
  
  // shuffle the row or column
  if (n == string("row"))
    random_shuffle(new_row_names.begin(), new_row_names.end(), myrandom);
  else if (n == string("col"))
    random_shuffle(new_col_names.begin(), new_col_names.end(), myrandom);
  else 
  {
    stringstream err;
    err << "ERROR: could not permute by " << n << ". options are row and col" << endl;
    error(err.str());
  }
  
  // permute the new data
  for (int i=0; i<nrows; i++)
  {
    string& rowname     = row_names[i];
    string& new_rowname = new_row_names[i];
    
    for (int j=0; j<ncols; j++)
    {
      string& colname     = col_names[j];
      string& new_colname = new_col_names[j];
      
      new_data[new_rowname][new_colname] = data[rowname][colname];
    }
  }
  
  // replace old data
  data = new_data;
  row_data.clear();
  col_data.clear();
  vectorize();
  
  // replace the node
  node->erase(name);
  write(*node, name, precision);
}
    
 
      
/*    I/O    */

template< typename T >
void DataTable<T>::read(ptree& pt, string node_name)
{
  ptree& table_node = pt.get_child(node_name);
 
  node = &pt;
  name = node_name;
  
  row_datatype = table_node.get<string>("<xmlattr>.row");
  col_datatype = table_node.get<string>("<xmlattr>.col");
  
  foreach_(ptree::value_type const& row, table_node)
  {
    if (row.first != "TableRow") continue;
    
    string row_name = row.second.get<string>(string("<xmlattr>." + row_datatype));
    if (!hasRowName(row_name)) row_names.push_back(row_name);
    //cerr << row_name << endl;
    const ptree & cols = row.second.get_child("<xmlattr>");
    foreach_(ptree::value_type const& col, cols)
    {
      if (col.first != row_datatype)
      {
        string col_name = col.first;
        if (!hasColName(col_name)) col_names.push_back(col_name);
        //cerr << col_name << endl;
        data[row_name][col_name] = row.second.get<T>(string("<xmlattr>." + col_name));
      }
    }   
  }
  fillNaN();
  vectorize();
}

template< typename T >
void DataTable<T>::set(vector<vector<T> > new_data, vector<string>& row_names, vector<string>& col_names, string row_datatype, string col_datatype)
{
  this->row_names = row_names;
  this->col_names = col_names;
  
  this->row_datatype = row_datatype;
  this->col_datatype = col_datatype;
  
  if (new_data.size() != row_names.size())
    error("row names must be same length as rows in data!");
  
      
  data.clear();
  
  int nrow = row_names.size();
  int ncol = col_names.size();
  for (int i=0; i<nrow; i++)
  {
    if (new_data[i].size() != col_names.size())
      error("col names must be same length as col in data!");
    const string& row_name = row_names[i];
    for (int j=0; j<ncol; j++)
    {
      const string& col_name = col_names[j];
      data[row_name][col_name] = new_data[i][j];
    }
  }
  vectorize();
}
  
  

template< typename T >
void DataTable<T>::write(ostream & os, string node_name, int precision)
{
  ptree pt;
  write(pt, node_name, precision);
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

template< typename T >
void DataTable<T>::write(ptree & pt, string node_name, int precision)
{
  int p = precision;
  int w = p + 7;
  int maxname = 0;
  int nrows = row_names.size();
  int ncols = col_names.size();
  
  // get the largest name for pretty printing
  for (int i=0; i<ncols; i++)
  {
    maxname = max(maxname, (int) col_names[i].size());
  }
  
  stringstream name;
  stringstream value;
  
  ptree & table_node = pt.add(node_name,"");
  table_node.put("<xmlattr>.row", row_datatype);
  table_node.put("<xmlattr>.col", col_datatype);
  
  for (int i = 0; i<nrows; i++)
  {
    ptree & row = table_node.add("TableRow","");
    row.put(string("<xmlattr>." + row_datatype), row_names[i]);
    for (int j=0; j<ncols; j++)
    {
      name.str("");
      value.str("");
      name << col_names[j];
      value << setw(w) << setprecision(p) << fixed << data[row_names[i]][col_names[j]];
      row.put("<xmlattr>."+name.str(),value.str());
    }
  }
}



// this function inverts the data table. Note this throws errors if IDs are numbers
/*
template< typename T >
void DataTable<T>::write(ptree & pt, string node_name)
{
  int w = 8;
  int precision = 4;
  int maxname = 0;
  int nrows = row_names.size();
  int ncols = col_names.size();
  
  // get the largest name for pretty printing
  for (int i=0; i<nrows; i++)
  {
    maxname = max(maxname, (int) row_names[i].size());
  }
  
  stringstream name;
  stringstream value;
  
  ptree & table_node = pt.add(node_name,"");
  table_node.put("<xmlattr>.col", row_datatype);
  table_node.put("<xmlattr>.row", col_datatype);
  
  for (int i = 0; i<ncols; i++)
  {
    ptree & row = table_node.add("TableRow","");
    row.put(string("<xmlattr>." + col_datatype), col_names[i]);
    for (int j=0; j<nrows; j++)
    {
      name.str("");
      value.str("");
      name << row_names[j];
      value << setw(6) << setprecision(4) << fixed << data[row_names[j]][col_names[i]];
      row.put("<xmlattr>."+name.str(),value.str());
    }
  }
}
*/  
  
template class DataTable<double>;
template class DataTable<int>;
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

