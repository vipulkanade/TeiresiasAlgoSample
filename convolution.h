#ifndef CONVOLUTION_H_
#define CONVOLUTION_H_

#include "pattern.h"

#include <map>
#include <vector>
#include <iostream>
#include <ostream>

using namespace std;

inline unsigned bracket_length_comp(const std::string motif1, const std::string motif2){
    size_t p1 = motif1.find('[');
    size_t p2 = motif2.find('[');
    //if no [ is encountered return true only if motif1 is longer
    if(p1 == string::npos && p2 == string::npos){
        if (motif1.size() > motif2.size()) return 1;
        else if (motif1.size() < motif2.size()) return 2;
        else return 0;
    }
    
    return 3;
}

struct prefix_wise_less
{
	bool operator()(const std::string &motif1, const std::string &motif2) const
	{
		return true;
	}
};

struct suffix_wise_less
{
	bool operator()(const std::string &motif1, const std::string &motif2) const
	{
		return true;
	}
};

struct Maximal_key
{
	size_t offset_list_size;
	long int diff_sum;
    unsigned first_sequence;
    unsigned last_sequence;
    long int first_diff_sum_sum;
    long int last_diff_sum_sum;
};


struct Max_ltOp
{
	bool operator()(const Maximal_key &lhs, const Maximal_key &rhs) const
	{
		if (lhs.offset_list_size < rhs.offset_list_size) return true;
        else if (lhs.offset_list_size > rhs.offset_list_size) return false;
		else if (lhs.diff_sum < rhs.diff_sum) return true;
        else if (lhs.diff_sum > rhs.diff_sum) return false;
        else if (lhs.first_sequence < rhs.first_sequence) return true;
        else if (lhs.first_sequence > rhs.first_sequence) return false;
        else if (lhs.last_sequence < rhs.last_sequence) return true;
        else if (lhs.last_sequence > rhs.last_sequence) return false;
        else if (lhs.first_diff_sum_sum < rhs.first_diff_sum_sum) return true;
        else if (lhs.first_diff_sum_sum > rhs.first_diff_sum_sum) return false;
        else if (lhs.last_diff_sum_sum < rhs.last_diff_sum_sum) return true;
        else if (lhs.last_diff_sum_sum > rhs.last_diff_sum_sum) return false;
        else return false;
	}
};

struct Soln_Op
{
	bool operator()(const Pattern &lhs, const Pattern &rhs) const
	{
		return lhs.Ls.ijs.size() > rhs.Ls.ijs.size();
	}
};

typedef std::pair<const std::string, Offset_list> Pattern_pair;
typedef const Pattern_pair *Ptr_pattern_pair;
typedef std::map<std::string, std::vector<Ptr_pattern_pair>, prefix_wise_less> Dir_p; // map from a prefix to a Ptr_pattern_pair
typedef std::map<std::string, std::vector<Ptr_pattern_pair>, suffix_wise_less> Dir_s; // map from a suffix to a Ptr_pattern_pair
typedef std::multimap<Maximal_key, Pattern, Max_ltOp> Maximal_map;

unsigned bracket_length_comp2(const std::string motif1, const std::string motif2);
bool dir_p_comp(std::string motif1, std::string motif2);
bool dir_s_comp(std::string motif1, std::string motif2);

// get the prefix or suffix of a pattern of size overlap using delim '.' for a pattern
std::string prefix(std::string motif, size_t overlap, char delim);
std::string suffix(std::string motif, size_t overlap, char delim);

// make the directory structures for eps organized by prefix and suffix
Dir_p make_dir_p(const std::vector<Elementary_patterns> &dynamic_eps, const size_t overlap_len);
Dir_s make_dir_s(const std::vector<Elementary_patterns> &dynamic_eps, const size_t overlap_len);

// FOR DEBUGGING ONLY: output directory structures
void out_dirs(std::ostream &out, const Dir_p &dir_p, const Dir_s &dir_s);

//Removes p from the dirs
void remove_entries(Pattern p, Dir_p& dir_p, Dir_s& dir_s, size_t overlap, char delim);

// left_convovlve q from Dir_s with t
Pattern left_convolve(Ptr_pattern_pair pQ, const Pattern &t, const std::string &prefix_w);

// right_convolve t with q from Dir_p
Pattern right_convolve(const Pattern &t, Ptr_pattern_pair pQ, const std::string &suffix_w);

//Count the chars without the brackets
unsigned non_bracket_diff_length(std::string first, std::string second);

int kth_character(int k,std::string s);

unsigned count_brackets(std::string s);

unsigned count_bracketed_length(std::string s);

// copy the ep from the Ptr_pattern_pair
inline Pattern copy_ep(Elementary_patterns::const_iterator itP)
{
	std::vector<IJ> ijs;
	Offset_list ol(ijs);
	Pattern ep(ol);
	ep.motif = itP->first;
	ep.Ls = itP->second;
	return ep;
}

// Creates the Maximal_key
Maximal_key make_max_key(const Pattern &pattern, const std::vector<Sequence> &seqs);

// modified trivial string search; modified KMP might be of value; returns true if pat_motif is a substr of max_motif
bool string_search(const Pattern &max_motif, const Pattern &pat_motif);


inline bool is_subsumed(const std::pair<Maximal_map::const_iterator, Maximal_map::const_iterator> ret, const Pattern &pattern,Maximal_map &maximal, Maximal_key key)
{
	for (Maximal_map::const_iterator it = ret.first; it != ret.second; ++it)
		if (string_search(it->second, pattern)) return true;
	return false;
}

// Creates the maximal key and determines whether or not to add the pattern into maximal
inline void add_pattern(Maximal_map &maximal, const Pattern &pattern, const std::vector<Sequence> &seqs)
{
	Maximal_key max_key = make_max_key(pattern, seqs);
	std::pair<Maximal_map::const_iterator, Maximal_map::const_iterator> ret = maximal.equal_range(max_key);
	if (ret.first == ret.second) maximal.insert(std::make_pair(max_key, pattern));
	else if (!is_subsumed(ret, pattern, maximal,max_key)) maximal.insert(std::make_pair(max_key, pattern));
}

// returns true if the key is not subsumed by something in the Maximal_map
inline bool is_maximal(Maximal_map &maximal, const Pattern &pattern, const std::vector<Sequence> &seqs)
{
	Maximal_key max_key = make_max_key(pattern, seqs);
	std::pair<Maximal_map::const_iterator, Maximal_map::const_iterator> ret = maximal.equal_range(max_key);
	if (ret.first == ret.second) return true;
	else if (!is_subsumed(ret, pattern, maximal,max_key)) return true;

    return false;
}

//Copies a stack of patterns to a vector of ones
std::vector<Pattern> clean_up_soln(const Maximal_map &maximal, std::vector<Sequence> seqs);
std::string clean_up_motif(Pattern p, std::vector<Sequence> seqs);
std::string clean_up_bracket(std::vector<IJ> offsets, long int motif_pos, std::vector<Sequence> seqs);

// outputs the solution vector according to config's specifications
void out_soln(std::vector<Pattern> &soln_vec, const Config &config, const std::vector<Sequence> &seqs);

Maximal_map remove_max_bracks(const Maximal_map maximal, unsigned n, vector<Sequence> seqs);
Maximal_map remove_big_k_patterns(const Maximal_map maximal, unsigned q);

#endif
