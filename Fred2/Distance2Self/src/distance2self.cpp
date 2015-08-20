#define BOOST_PYTHON_STATIC_LIB
#include <boost/python.hpp>
#include "distance2self.hpp"

using namespace std;
using namespace boost::python;


Distance2Self::Distance2Self(string f_m, boost::python::list selfs, int peptide_length):
    m_(f_m), ta_()
{
    //Matrix m("/abi-projects/dist2self/matrices/BLOSUM45_distance_normal.dat"); cout << "Initializing trie. " << endl; Trie t;
    Trie t;
    Matrix::IndexSequence indices;
    for (int i = 0; i < len(selfs); ++i)
    {
        string s = boost::python::extract<string>(selfs[i]);
        m_.translate(s, indices);
        t.add(indices);
    }
    ta_ = TrieArray(t, peptide_length);
}

d2s_result Distance2Self::getDistance2Self(boost::python::list sequences)
{
    map<string, map<string,double> > peptides;
    for (int i = 0; i < len(sequences); ++i)
    {
        //http://stackoverflow.com/questions/17659572/boostpython-passing-reference-of-pythonlist
        string s = boost::python::extract<string>(sequences[i]);
        peptides.insert(make_pair(s,map<string,double>()));
    }

    for( map<string,map<string,double> >::iterator iter = peptides.begin(); iter != peptides.end(); ++iter )
    {
        vector<size_t> p, seq_i;
        m_.translate(iter->first, seq_i); //translate peptide sequence to matrix indices. seq contains the translated peptide sequence.

        multiset<pair<double,string> > dist_min, dt;
        double dist = 0.0;
        pair<size_t, size_t> n (0,0); //Node (0,0) : start at top of the trie
        dist_min = DFS_BnB_x_pair(ta_,n,dist,m_,p,seq_i,dt);
        for (multiset<pair<double,string> >::iterator dists=dist_min.begin() ; dists != dist_min.end(); dists++ )
        {
            peptides[iter->first][dists->second] = dists->first;
        }
    }
    return peptides;
}

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
    class_<sequence_list>("sequence_list")
        .def(vector_indexing_suite<sequence_list>() );

    class_<sequence_map>("sequence_map")
        .def(map_indexing_suite<sequence_map>() );

    class_<d2s_result>("d2s_result")
        .def(map_indexing_suite<d2s_result>() );

    class_<Distance2Self>("Distance2Self",init<std::string, boost::python::list, int>())
            .def("setTrieArray", &Distance2Self::setTrieArray)
            .def("getDistance2Self", &Distance2Self::getDistance2Self)
            .def("setTrieArray", &Distance2Self::setTrieArray)
            .def("setMatrix", &Distance2Self::setMatrix);
}
