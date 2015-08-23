#ifndef DISTANCE2SELF_H
#define DISTANCE2SELF_H

#include <boost/python.hpp>
#include <vector>
#include <string>
#include "TrieArray.hpp"

using namespace std;

class Distance2Self
{
  private:
    Distance2Self();
  protected:
    Matrix m_;
    TrieArray ta_;

  public:
    Distance2Self(string f_m, boost::python::list selfs);
    boost::python::dict getDistance2Self(boost::python::list sequences, int max_result);
    void setTrieArray(string f_ta);
    void setMatrix(string f_m);
//    Distance2Self();
};

#endif // DISTANCE2SELF_H
