#ifndef PATTERN_H_
#define PATTERN_H_

#include "config.h"
#include "seqs.h"
#include "alphabet.h"

#include <string>
#include <map>
#include <vector>

using namespace std;


struct IJ
{
	unsigned i;
	unsigned j;
};

struct Offset_list
{
	Offset_list(const std::vector<IJ>& aijs);
    std::vector<IJ> ijs;
};

struct Pattern
{
	Pattern(const Offset_list& ol);
	std::string motif;
	Offset_list Ls;
};

struct alphabetical
{
	bool operator()(const std::string &motif1, const std::string &motif2) const
	{
		size_t len = 0;
		motif1.size() <= motif2.size() ? len = motif1.size() : len = motif2.size();
		for (size_t i = 0; i < len; ++i)
		{ // density preferred left to right
            if ((motif1.at(i) == '.') && (motif2.at(i) == '[')) return false; //dot is prefix-wise less than bracket
            else if ((motif1.at(i) == '[') && (motif2.at(i) == '.')) return true;
			else if ((motif1.at(i) != '.') && (motif2.at(i) == '.')) return true; // motif1 is prefix-wise less
			else if ((motif1.at(i) == '.') && (motif2.at(i) != '.')) return false;
            else if ((motif1.at(i) != '[') && (motif2.at(i) == '[')) return true;
            else if ((motif1.at(i) == '[') && (motif2.at(i) != '[')) return false;
		}
		if (motif1.size() >= len) return true; // motif1 is longer
		else if (motif1.size() < len) return false;
		else if (motif1.compare(motif2) < 0) return true; // ties resolved alphabetically
		else return false;
	}
};

typedef std::map<std::string, Offset_list, alphabetical> Elementary_patterns;

// makes the starting bit_mask and returns it
std::string initialize_bit_mask(const unsigned w, const unsigned l);
std::string initialize_int_mask(const unsigned w, const unsigned l, const unsigned max_b);


// trim trailing 0s from bit_mask
std::string trim0s(const std::string& bit_mask);

//
std::string make_and_insert(charlist substs, Elementary_patterns& eps, const string& current_bit_mask, const vector<Sequence>& seqs,const unsigned i,const unsigned j);

// checks if the eps are supported
Elementary_patterns supported(Elementary_patterns& eps, const Config& config);

// outputs the eps from the scanning phase
void out_scan(const std::vector<Elementary_patterns>& dynamic_eps, const Config& config, const std::vector<Sequence>& seqs);

#endif
