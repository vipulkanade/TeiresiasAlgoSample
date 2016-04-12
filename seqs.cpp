#include "seqs.h"

#include <fstream>
#include <iostream>

using namespace std;

inline void skip_empty_lines(ifstream *pcinfs)
{
	while(isspace(pcinfs->peek()))
		pcinfs->get();
}


inline void trim(string *pstr)
{
	while(isspace((*pstr)[pstr->size() - 1]))
		pstr->erase(pstr->end() - 1);
}

vector<Sequence> parse_seqs(const string& file_name)
{
	ifstream cinfs(file_name.c_str());
	vector<Sequence> seqs;
    vector<Sequence>::iterator it;
    Sequence current_seq;
    bool header=false;
	if (cinfs.is_open())
	{
		const char DELIM_START = '>';
		const char DELIM_END = ' ';
		if (!cinfs.good()) cerr << "File stream is not good." << endl;
		while (cinfs.good())
		{
			skip_empty_lines(&cinfs);
            if (cinfs.peek() == DELIM_START){
                cinfs.get();
                getline(cinfs, current_seq.header);
                header = true;
            }
            else{
                getline(cinfs, current_seq.str);
                trim(&current_seq.str);
                if(seqs.end()==seqs.begin()){
                    if (header==false) current_seq.header="No header provided";
                    header=true;
                }
                if(header==true){
                    seqs.push_back(current_seq);
                }
                else{
                    it = seqs.end();
                    it--;
                    it->str += current_seq.str;
                }
                current_seq.header.erase(0,string::npos);
                current_seq.str.erase(0,string::npos);
                header = false;
            }
		}
	} else cerr << "Can't open input file." << endl;
	return seqs;
}

void out_seqs(ostream& out, const vector<Sequence>& seqs)
{
	for (vector<Sequence>::const_iterator iter = seqs.begin(); iter != seqs.end(); ++iter)
	{
		out << "Header: " << iter->header << endl;
		out << "String: " << iter->str << endl;
	}
}

size_t smallest_seq(const vector<Sequence>& seqs)
{
	size_t length = 0;
	for (vector<Sequence>::const_iterator iter = seqs.begin(); iter != seqs.end(); ++iter)
	{
		if ((iter->str.size() < length) || (length == 0))
			length = iter->str.size();
	}
	return length;
}

Sequence_line transform_seq (vector<Sequence> seqs)
{
    Sequence_line seqsl;
    for (vector<Sequence>::const_iterator iter = seqs.begin(); iter != seqs.end(); ++iter){
        seqsl.seq_len.push_back(iter->str.size());
        seqsl.sequences += iter->str;
    }
    return seqsl;
}
