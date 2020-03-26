/*********************************************************************************
*                                                                                *
*     utils.h                                                                    *
*                                                                                *
*     Contains some general functions useful for a number of different classes   *
*                                                                                *
*********************************************************************************/

#include <vector>
#include <string>
#include <iostream>
#include <limits>
#include <boost/shared_ptr.hpp>

using namespace std;

// convert ACGTN to 01234;
vector<int> string2int(const string& s);

string int2string(const vector<int> & s);

vector<char> int2char(const vector<int>& s);

vector<int> char2int(const vector<char>& s);

// error message

void error(const string & error_message);

void warning(const string & warn_message);

double correlation(vector<double>& x, vector<double>& y);

//double sse(vector<double*>& x, vector<double*>& y);
