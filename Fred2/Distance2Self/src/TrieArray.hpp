#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <tr1/cstdint>
#include "Trie.hpp"
#include "Matrix.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

using namespace std;

class TrieArray
{
  public:
    typedef size_t PackedIndex;
    typedef vector<PackedIndex> LevelArray;
    typedef vector<LevelArray> PackedTrie;

    PackedTrie data;
	TrieArray();
    TrieArray(string ta);
    TrieArray(const Trie& t);
 	void add(string s);
	
	template<class Archive>
		void save(Archive & ar, const unsigned int version) const;
		
	template<class Archive>
		void load(Archive & ar, const unsigned int version);

    BOOST_SERIALIZATION_SPLIT_MEMBER();

  private:
	friend class boost::serialization::access;

};
BOOST_CLASS_VERSION(TrieArray, 1);


TrieArray::TrieArray()
    : data ()
{
}

TrieArray::TrieArray(string ta)
    : data ()
{
    std::ifstream ifs(ta.c_str());
    boost::archive::text_iarchive ia(ifs);
    load(ia,1);
}


TrieArray::TrieArray(const Trie& t)
{
    vector<Trie::TrieGraph::vertex_descriptor> nodes;
    nodes.push_back(t.root);

    data.push_back(LevelArray(1, 0));
    while (!nodes.empty())
    {
        size_t number_of_nodes = 0;
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            data.back()[i] |= number_of_nodes;
            number_of_nodes += out_degree(nodes[i], t.g);
        }
        data.push_back(LevelArray());
        data.back().resize(number_of_nodes);
        //cout << "Layer " << data.size() << " has " << number_of_nodes << " nodes." << endl;
        size_t j = 0;
        vector<Trie::TrieGraph::vertex_descriptor> new_nodes(number_of_nodes);
        for (size_t i = 0; i < nodes.size(); ++i)
        {
            Trie::TrieGraph::adjacency_iterator ai, a_end;
            for (tie(ai, a_end) = adjacent_vertices(nodes[i], t.g); ai != a_end; ++ai, ++j)
				{
					data.back()[j] = get(vertex_index_t(), t.g, *ai) << 56;
					new_nodes[j] = *ai;
				}
        }
        nodes = new_nodes;
    }
}
  
void TrieArray::add(string s)
{
    //cout << "Adding: " << s << endl;
    //TODO
}

template<class Archive>
void TrieArray::save(Archive & ar, const unsigned int version) const
{
    for (size_t i=0; i<data.size(); ++i)
  	{
  						cout << data.size() << endl;
						ar & data[i];
  	} 	
}

template<class Archive>
void TrieArray::load(Archive & ar, const unsigned int version)
{
    for (size_t i=0; i<11; ++i)
    {
        LevelArray v;
        ar & v;
        data.push_back(v);
    }
}


typedef pair<size_t, size_t> Node;
typedef vector<Node> Children;
typedef vector<size_t> Peptide;

static Children getChildren(const TrieArray& ta, Node n)
{
    Children children;
    if (n.first == ta.data.size()-1) //leaf
    {
        return children;
    }
    else
    {
        size_t v_e;

        size_t v = ta.data[n.first][n.second]; //Inhalt aktueller Knoten
        size_t v_s = (v & 0x00ffffffffffffff); //Zeiger auf erstes Kind bzw. Index des ersten Kindes in naechster Ebene.

        if (n.second >= ta.data[n.first].size()-1) //last sibling
        {
            v_e = ta.data[n.first + 1].size();
        }
        else
        {
            size_t v_next = ta.data[n.first][n.second+1]; //Inhalt naechster Geschwisterknoten, braucht man fuer den Zeiger auf Ende der Kinder vom aktuellen Knoten
            v_e = (v_next & 0x00ffffffffffffff); //Zeiger auf letztes Kind bzw. Index des letzten Kindes in naechster Ebene.
        }

        //if (v_s > v_e)
        //{
        //    cout << v << " " << v_s << " " << v_e << " n1: " << n.first << " n2: " << n.second << " " << ta.data[n.first].size() << endl;
        //}

        for( size_t i = v_s; i < v_e; i++ )
        {
            size_t level = n.first + 1;
            size_t sibling = i;
            Node c (level,sibling);
            children.push_back(c);
        }
        return children;
    }
}


static multiset<pair<double, string> > DFS_BnB_x_pair(const TrieArray& ta, Node n, double& dist, Matrix& m, Peptide& p, Peptide& s, multiset<pair<double, string> >& dist_min, size_t x=10)
/*
Returns the x clostest peptides
Returns pairs of distance and peptide
Parameters:
TrieArray: The search trie
dist: distance of search peptide s to
m: Distance matrix //we have to see what we will do with similarity matrices...
p: path to the current node in the trie. keeps track of the words stored in the nodes.
s: search peptide
dist_min: minimal distance foud so far. Required for break in branch and bound.
*/
{
    if (dist_min.size() >= x) //If not enough distances found so far
    {
        if (dist > (*--dist_min.end()).first)//check if the bound criteria is fullfilled (distance is already greater or equal to maximal distance selected so far)
        {
            //cout << "break: " << *--dist_min.end() << " (" << dist << ") " <<  "-----";
            return dist_min;
        }
    }

    Children c = getChildren(ta, n);
    if (c.size() == 0) //at leaf
    {
        string seq;
        for (size_t pos=0; pos<p.size();++pos)
        {
            seq += m.indexToChar(p[pos]);
        }
        if (dist_min.size() < x)
        {
            pair<double, string> temp_pair(dist,seq);
            dist_min.insert(temp_pair);
            //cout << dist << "\t" << seq << endl;
        }
        else if (dist < (*--dist_min.end()).first)//peptide with new minimal distance found.
        {
            pair<double, string> temp_pair(dist,seq);
            //cout << seq << endl;
            dist_min.insert(temp_pair);
            //cout << dist << "\t" << seq << endl;
            dist_min.erase(--dist_min.end()); //remove largest distance from list
            //cout << --dist_min.end() << endl;
        }

        return dist_min;
    }

    for (size_t i=0; i<c.size(); ++i)
    {
        size_t a = (ta.data[c[i].first][c[i].second] >> 56);
        size_t b = s[(c[i].first-1)];
        p.push_back(a);
        double d  = m.data_[a*20 + b];
        if (abs(d)<0.0000001)
            {d = 0.0;}
        dist += d;
        multiset<pair<double,string> > dist_min_tmp;
        dist_min_tmp = DFS_BnB_x_pair(ta,c[i],dist,m,p,s, dist_min, x);
        dist_min = dist_min_tmp;
        p.pop_back();
        dist -= d;
        if (abs(dist)<0.0000001)
            {dist = 0.0;}
    }
    return dist_min;
}

