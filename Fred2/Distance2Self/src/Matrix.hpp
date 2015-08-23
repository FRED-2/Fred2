#include <fstream>

using namespace std;

class Matrix
{
  public:
	typedef vector<size_t> IndexSequence;

  protected:
    vector<size_t> matrix_;
    size_t alphabet_size_;
    string alphabet_;
    map<char, size_t> char_to_index_;
    map<size_t, char> index_to_char_;

  public:
    Matrix(const string& filename);
    vector<double> data_;

	size_t charToIndex(char c) { return char_to_index_[c]; }
	char indexToChar(size_t i) { return index_to_char_[i]; }
	
	void translate(const string& s, IndexSequence& sequence);
  protected:
	void setupMaps_();
  private:
    Matrix();
};

Matrix::Matrix(const string& filename)
	: alphabet_(""), data_(20 * 20)
{
    //TODO: check file exists
    ifstream is(filename.c_str());
	data_.resize(20 * 20);
	for (size_t i = 0; i < 20; ++i)
	{
		char c;
		is >> c;
		alphabet_.push_back(c);
		double score;
		for (size_t j = 0; j < 20; ++j)
		{
			is >> score;
			data_[i * 20 + j] =  score;
		}
	}
	setupMaps_();
}

void Matrix::setupMaps_()
{
	for (size_t i = 0; i < alphabet_.size(); ++i)
	{
		index_to_char_.insert(make_pair(i, alphabet_[i]));
		char_to_index_.insert(make_pair(alphabet_[i], i));
	}
}

// Matrices are stored with integers as indices, not with letters
void Matrix::translate(const string& s, Matrix::IndexSequence& sequence)
{
    sequence.resize(s.size());
    for (size_t i = 0; i < s.size(); ++i)
    {
        sequence[i] = char_to_index_[s[i]];
	}
}

