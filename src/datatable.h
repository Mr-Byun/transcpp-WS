/*********************************************************************************
*                                                                                *
*     datatable.h                                                                *
*                                                                                *
*     This contains a simple data-table class. It holds a 2D table and provides  *
*     methods for access and modification. Each row and column has names used    *
*     for access. This also contains methods for reading and writing tables to   *
*     fixed width tables or xml tables. This uses a map implementation with      *
*     string access, so for efficiency inner loop calculations should act on     *
*     pointers to the data contained here rather than calling the get functions  *
*                                                                                *
*********************************************************************************/

#ifndef DATATABLE_H
#define DATATABLE_H

#include "mode.h"

#include <vector>
#include <string>
#include <map>
#include <cstdlib>
#include <boost/unordered_map.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <limits>
#include <iostream>

using namespace std;
using boost::property_tree::ptree;

/* the data in every row and column must be the same data type in this simple
implementation, which should be sufficient for everything we need it for */

/* Note that map access by string is slow! Dont use this in the inner loop! */

template< typename T> 
class DataTable
{
private: 
  mode_ptr mode;
  ptree* node;
  string name;
  
  // keep track of which variable is in rows and which is in columns
  string row_datatype;
  string col_datatype;
  
  // row and column names
  vector<string> row_names;
  vector<string> col_names;
  
  // the primary data, using data[rowname][colname]
  map<string, map<string, T> > data;
  
  // vectorized rows and columns for quick access
  map<string, vector<T*> > row_data;
  map<string, vector<T*> > col_data;
  
public:
  // Constructors
  DataTable();
  DataTable(ptree& pt, string name);
  DataTable(string fname, string name, string section);
  
  // Getters
  T& getDataPoint(string type1, string id1, string type2, string id2); 
  vector<T*>& getRow(string n);
  vector<T*>& getCol(string n);
  const string& getRowType() { return row_datatype; }
  const string& getColType() { return col_datatype; }
  int getN(string);
  int getNrows() { return row_names.size(); }
  int getNcols() { return col_names.size(); }
  vector<string>& getRowNames() { return row_names; }
  vector<string>& getColNames() { return col_names; }
  vector<string>& getNames(string type);
  
  // Setters
  void setMode(mode_ptr mode) { this->mode = mode; }
  void set(vector<vector<T> > data, vector<string>& row_names, vector<string>& col_names, string row_datatype, string col_datatype);
  
  // methods
  bool hasRowName(string& n);
  bool hasColName(string& n);
  void fillNaN();
  void vectorize();
  void permute(string& n, int precision); // n is "row" or "col" for which to permute 

  // Output
  void read(ptree& pt, string name);
  void write(ostream & os, string node_name, int precision);
  void write(ptree & pt, string node_name, int precision);
};

typedef boost::shared_ptr<DataTable<double> > table_ptr;


























#endif
