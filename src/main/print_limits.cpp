#include <iostream>
#include <iomanip>
#include <limits>
#include <cstdlib>

using namespace std;
int main(int argc, char* argv[])
{
  cerr << setw(20) << "type"        << setw(20) << "max"                              << endl;
  cerr << setw(20) << "float"       << setw(20) << numeric_limits<float>::max()       << endl;
  cerr << setw(20) << "double"      << setw(20) << numeric_limits<double>::max()      << endl;
  cerr << setw(20) << "long double" << setw(20) << numeric_limits<long double>::max() << endl;
  return 0;
}
  