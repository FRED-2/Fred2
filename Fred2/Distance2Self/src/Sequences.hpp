class Sequences
{
	public:
	typedef vector<string> SequenceVector;
	

	protected:
	SequenceVector data_;
	SequenceVector desc_;

	public:
	string& operator [] (size_t index) { return data_[index]; }
	const string& operator [] (size_t index) const { return data_[index]; }
	
	string& description(size_t index) { return desc_[index]; }
	const string& description(size_t index) const { return desc_[index]; }

	size_t size() const { return data_.size(); }

	void clear() { data_.resize(0); }

	Sequences(size_t size = 0) : data_(size), desc_(size) {}

	Sequences(const string& filename)
	{

		ifstream is(filename.c_str());
		if (not is)
			throw "Cannot open File!";

    string line;
    bool done = false;

    while	(not done and getline(is, line)) 
		{
      // quit reading after a blank line
      if (line.size() == 0) break;

      // compare if expectations are met...
      if (line[0] != '>')
        throw "FASTA sequence doesn't start with '>'";

      // Parse the header
      string desc = line;

      // Parse the sequence
      string s;
      while (is) 
			{
        char c = is.get();

        // bail on EOF
        if (not is) { done = true; break;}
        is.putback(c);

        // finish this sequence before beginning of next sequence
        if (c == '>') break;

        // read the next line of letters
        getline(is, line);

        // quit reading after a blank line
        if (not line.size()) { done = true; break;}

        // add the letters in
        s += line;
			}

      data_.push_back(s);
			desc_.push_back(desc);
		}
	}
	
	vector<string>::iterator getDataBegin()
		{
			return data_.begin();
		}
		
	
	vector<string>::iterator getDataEnd()
		{
			return data_.end();
		}
		
	void push_back(const string& s) 
	{
		data_.push_back(s);
		desc_.push_back("");
	}

	void push_back(const string& s, const string& d) 
	{
		data_.push_back(s);
		desc_.push_back(d);
	}
	
	
};
