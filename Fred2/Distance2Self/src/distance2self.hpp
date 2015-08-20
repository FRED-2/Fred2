#ifndef DISTANCE2SELF_H
#define DISTANCE2SELF_H

#include <vector>
#include <string>
#include "TrieArray.hpp"

using namespace std;

typedef vector<string> sequence_list;
typedef map<string, double> sequence_map;
typedef map<string, sequence_map > d2s_result;

class Distance2Self
{
  private:
    Distance2Self();
  protected:
    Matrix m_;
    TrieArray ta_;

  public:
    Distance2Self(string f_m, boost::python::list selfs, int peptide_length);
    d2s_result getDistance2Self(boost::python::list sequences);
    void setTrieArray(string f_ta);
    void setMatrix(string f_m);
//    Distance2Self();
};

#endif // DISTANCE2SELF_H
