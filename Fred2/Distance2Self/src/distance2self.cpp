#define BOOST_PYTHON_STATIC_LIB
#include "distance2self.hpp"

using namespace std;
using namespace boost::python;

Distance2Self::Distance2Self(string f_m, boost::python::list selfs):
    m_(f_m), ta_()
{
    //Matrix m("/abi-projects/dist2self/matrices/BLOSUM45_distance_normal.dat"); cout << "Initializing trie. " << endl; Trie t;
    Trie t;
    Matrix::IndexSequence indices;
    for (int i = 0; i < len(selfs); ++i)
    {
        //TODO check for equi-legth?!
        string s = boost::python::extract<string>(selfs[i]);
        m_.translate(s, indices);
        t.add(indices);
    }
    ta_ = TrieArray(t);
}

boost::python::dict Distance2Self::getDistance2Self(boost::python::list sequences, int max_result)
{
    dict peptides;
    for (int i = 0; i < len(sequences); ++i)
    {
        //http://stackoverflow.com/questions/17659572/boostpython-passing-reference-of-pythonlist
        string seq = boost::python::extract<string>(sequences[i]);

        vector<size_t> p, seq_i;
        m_.translate(seq, seq_i); //translate peptide sequence to matrix indices. seq contains the translated peptide sequence.

        multiset<pair<double,string> > dist_min, dt;
        double dist = 0.0;
        pair<size_t, size_t> n (0,0); //Node (0,0) : start at top of the trie
        dist_min = DFS_BnB_x_pair(ta_,n,dist,m_,p,seq_i,dt, max_result);
        dict dists;
        for (multiset<pair<double,string> >::iterator dit=dist_min.begin() ; dit != dist_min.end(); ++dit)
        {
            dists[dit->second] = dit->first;
        }
        peptides[seq]=dists;
    }
    return peptides;
}

//~ template <class K, class V>
//~ dict map2BoostPythonDict(std::map<K, V> map) {
    //~ typename map<K, V>::iterator it;
    //~ dict bpd;
    //~ for (it = map.begin(); it != map.end(); ++it) {
        //~ bpd[it->first] = it->second;
    //~ }
    //~ return bpd;
//~ }

void Distance2Self::setTrieArray(string f_ta)
{
    ta_ = TrieArray(f_ta);
}

void Distance2Self::setMatrix(string f_m)
{
    m_ = Matrix(f_m);
}

#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
BOOST_PYTHON_MODULE(d2s)
{
    class_<Distance2Self>("Distance2Self",init<std::string, boost::python::list>())
            .def("setTrieArray", &Distance2Self::setTrieArray)
            .def("getDistance2Self", &Distance2Self::getDistance2Self)
            .def("setTrieArray", &Distance2Self::setTrieArray)
            .def("setMatrix", &Distance2Self::setMatrix);
}
