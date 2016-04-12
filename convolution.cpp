#include "convolution.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>

using namespace std;


/* Input: A pattern, the prefix size (def: l-1) and one non literal char (def: '.')  */
/* Output: The prefix                                                                */
string prefix(string pattern, size_t overlap, char delim)
{
	string prefix = "";
	unsigned count_l = 0;
    size_t i;
	for (i = 0; count_l < overlap; ++i)
	{
        if(pattern.at(i) == '['){
            while(pattern.at(i) != ']'){
                prefix.push_back(pattern.at(i));
                i++;
            }
            prefix.push_back(pattern.at(i));
            ++count_l;
        }
		else if (pattern.at(i) != delim) {
			prefix.push_back(pattern.at(i));
			++count_l;
		}
        else{
            prefix.push_back(pattern.at(i));
        }
    }
	return prefix;
}

/* Input: A pattern, the suffix size (def: l-1) and one non literal char (def: '.')  */
/* Output: The suffix                                                                */
string suffix(string pattern, size_t overlap, char delim)
{
	string suffix = "";
	unsigned count_l = 0;
    size_t i;
	for (i = pattern.size() - 1; count_l < overlap; i--)
	{
        if(pattern.at(i) == ']'){
            while(pattern.at(i) != '['){
                suffix.push_back(pattern.at(i));
                i--;
            }
            suffix.push_back(pattern.at(i));
            ++count_l;
        }
        else if (pattern.at(i) != delim)
        {
            suffix.push_back(pattern.at(i));
            ++count_l;
        }
		else suffix.push_back(pattern.at(i));
	}
	reverse(suffix.begin(), suffix.end());
	return suffix;
}

unsigned bracket_length_comp2(const std::string motif1, const std::string motif2){
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

bool dir_p_comp(std::string motif1, std::string motif2){
    size_t len = 0;
    if(motif1.size() <= motif2.size())
        len = motif1.size();
    else
        len = motif2.size();
    size_t i;
    for (i = 0; i < len; ++i)
    { // density preferred left to right
        if ((motif1.at(i) == '.') && (motif2.at(i) == '[')) return false; //dot is prefix-wise less than bracket
        else if ((motif1.at(i) == '[') && (motif2.at(i) == '.')) return true;
        else if ((motif1.at(i) != '.') && (motif2.at(i) == '.')) return true; // motif1 is prefix-wise less
        else if ((motif1.at(i) == '.') && (motif2.at(i) != '.')) return false;
        else if ((motif1.at(i) != '[') && (motif2.at(i) == '[')) return true;
        else if ((motif1.at(i) == '[') && (motif2.at(i) != '[')) return false;
    }
    unsigned flag = bracket_length_comp2(motif1,motif2);
    if (flag == 1) return true; // motif1 is longer
    else if (flag == 2) return false;
    else if (motif1.compare(motif2) < 0) return true; // ties resolved alphabetically
    else return false;
}

bool dir_s_comp(std::string motif1, std::string motif2){
    size_t len = 0, len1 = motif1.size() - 1, len2 = motif2.size() - 1;
    if(motif1.size() < motif2.size())
        len = motif1.size();
    else
        len = motif2.size();
    size_t i;
    for (i = 0; i < len; ++i)
    { // density preferred right to left
        if ((motif1.at(len1 - i) == '.') && (motif2.at(len2 - i) == ']')) return false;
        else if ((motif1.at(len1 - i) == ']') && (motif2.at(len2 - i) == '.')) return true;
        else if ((motif1.at(len1 - i) != '.') && (motif2.at(len2 - i) == '.')) return true; // motif1 is suffix-wise less
        else if ((motif1.at(len1 - i) == '.') && (motif2.at(len2 - i) != '.')) return false;
        else if ((motif1.at(len1 - i) == ']') && (motif2.at(len2 - i) != ']')) return false;
        else if ((motif1.at(len1 - i) != ']') && (motif2.at(len2 - i) == ']')) return true;
    }
    if (len < motif1.size()) return true; // motif1 is longer
    else if (motif1.size() < len) return false;
    else if (motif1.compare(motif2) < 0) return true; // ties resolved alphabetically
    else return false;
}

void  extract_p(const Elementary_patterns& eps, const size_t overlap_len, Dir_p& dir_p)
{
    Elementary_patterns::const_iterator iter;
	for (iter = eps.begin(); iter != eps.end(); ++iter)
	{
		string ep_prefix = prefix(iter->first, overlap_len, '.');
		Ptr_pattern_pair pep = &(*iter);

		Dir_p::iterator itDirP = dir_p.begin();
        while(itDirP != dir_p.end()){
            if(itDirP->first == ep_prefix)
                break;
            itDirP++;
        }
        if (itDirP != dir_p.end()){
            vector<Ptr_pattern_pair>::iterator itVector;
            bool flag = false;
            for(itVector = itDirP->second.begin(); itVector != itDirP->second.end(); itVector++){
                if(dir_p_comp(pep->first, (*itVector)->first)){
                    itDirP->second.insert(itVector,pep);
                    flag = true;
                    break;
                }
            }
            if(!flag)
                itDirP->second.push_back(pep);
        }
		else
		{
			vector<Ptr_pattern_pair> pep_vec;
			pep_vec.push_back(pep);
			dir_p.insert(make_pair(ep_prefix, pep_vec));
		}
	}
}

inline void out_patterns(ostream& out, vector<Ptr_pattern_pair> pep_vec)
{
    vector<Ptr_pattern_pair>::const_iterator itP;
	for (itP = pep_vec.begin(); itP != pep_vec.end(); ++itP)
		out << (*itP)->first << " ";
}

inline void out_dir_p(ostream& out, const Dir_p& dir_p)
{
    out << "dir_p:" << endl;
    Dir_p::const_iterator iter;
	for (iter = dir_p.begin(); iter != dir_p.end(); ++iter)
	{
		out << iter->first << "->{ "; // prefix
		out_patterns(out, iter->second);
		out << "}" << endl;
	}
}

inline void out_dir_s(ostream& out, const Dir_s& dir_s)
{
    out << "dir_s:" << endl;
    Dir_s::const_iterator iter;
	for (iter = dir_s.begin(); iter != dir_s.end(); ++iter)
	{
		out << iter->first << "->{ "; // suffix
		out_patterns(out, iter->second);
		out << "}" << endl;
	}
}

Dir_p make_dir_p(const vector<Elementary_patterns>& dynamic_eps, const size_t overlap_len)
{
	Dir_p dir_p;
    vector<Elementary_patterns>::const_iterator iter;
	for (iter = dynamic_eps.begin(); iter != dynamic_eps.end(); ++iter)
		extract_p(*iter, overlap_len, dir_p);
    
    //out_dir_p(cout,dir_p);
	return dir_p;
}

void extract_s(const Elementary_patterns& eps, const size_t overlap_len, Dir_s& dir_s)
{
    Elementary_patterns::const_iterator iter;
	for (iter = eps.begin(); iter != eps.end(); ++iter)
	{
		string ep_suffix = suffix(iter->first, overlap_len, '.');
		Ptr_pattern_pair pep = &(*iter);
		Dir_s::iterator itDirS = dir_s.begin();
        while(itDirS != dir_s.end()){
            if(itDirS->first == ep_suffix)
                break;
            itDirS++;
        }
        if (itDirS != dir_s.end()){    
            vector<Ptr_pattern_pair>::iterator itVector;
            bool flag = false;
            for(itVector = itDirS->second.begin(); itVector != itDirS->second.end(); itVector++){
                if(dir_s_comp(pep->first, (*itVector)->first)){
                    itDirS->second.insert(itVector,pep);
                    flag = true;
                    break;
                }
            }
            if(!flag)
                itDirS->second.push_back(pep);
        }
		else
		{
			vector<Ptr_pattern_pair> pep_vec;
			pep_vec.push_back(pep);
			dir_s.insert(make_pair(ep_suffix, pep_vec));
		}
	}
}

Dir_s make_dir_s(const std::vector<Elementary_patterns>& dynamic_eps, const size_t overlap_len)
{
	Dir_s dir_s;
    vector<Elementary_patterns>::const_iterator iter;
	for (iter = dynamic_eps.begin(); iter != dynamic_eps.end(); ++iter)
		extract_s(*iter, overlap_len, dir_s);
    
    //out_dir_s(cout,dir_s);
	return dir_s;
}

void out_dirs(std::ostream& out, const Dir_p& dir_p, const Dir_s& dir_s)
{
	out_dir_p(out, dir_p);
	out << endl;
	out_dir_s(out, dir_s);
};

void remove_entries(Pattern p, Dir_p& dir_p, Dir_s& dir_s, size_t overlap, char delim){
    //make prefix
    string s = prefix(p.motif, overlap, delim);
    //Find and erase entry
    Dir_p::iterator itDir = dir_p.begin();
    while(itDir != dir_p.end()){
        if(itDir->first == s)
            break;
        itDir++;
    }
    vector<Ptr_pattern_pair>::iterator it = itDir->second.begin();
    while(it != itDir->second.end()){
        if((*it)->first == p.motif)
            break;
        it++;
    }
    itDir->second.erase(it);
    if (itDir->second.empty())
        dir_p.erase(itDir);

    //same for suffix
    s = suffix(p.motif, overlap, delim);
    itDir = dir_s.begin();
    while(itDir != dir_s.end()){
        if(itDir->first == s)
            break;
        itDir++;
    }
    it = itDir->second.begin();
    while(it != itDir->second.end()){
        if((*it)->first == p.motif)
            break;
        it++;
    }
    itDir->second.erase(it);
    if (itDir->second.empty())
        dir_s.erase(itDir);
}

            
            /*  CONVOLUTION  */


void match_ij(vector<IJ>::const_iterator itPIJ, vector<IJ>::const_iterator itSIJ, unsigned delta_j, Offset_list *pol)
{
	if (itPIJ->i == itSIJ->i)
	{
		if (itPIJ->j + delta_j == itSIJ->j)
		{
			IJ ij = {itPIJ->i, itPIJ->j};
			pol->ijs.push_back(ij);
		}
	}
}


// Just a trivial linear search
void match_ijs_lhs(Ptr_pattern_pair pQ, const Pattern& t, unsigned delta_j, Offset_list *pol)
{
    //iterator for the offlist of pattern pQ
    vector<IJ>::const_iterator itPIJ;
    vector<IJ>::const_iterator itSIJ;
	for (itPIJ = pQ->second.ijs.begin(); itPIJ != pQ->second.ijs.end(); ++itPIJ){
        //iterator for the offlist of pattern t
		for (itSIJ = t.Ls.ijs.begin(); itSIJ != t.Ls.ijs.end(); ++itSIJ){
			match_ij(itPIJ, itSIJ, delta_j, pol);
        }
	}
}

// Just a trivial linear search
void match_ijs_rhs(const Pattern& t, Ptr_pattern_pair pQ, unsigned delta_j, Offset_list *pol)
{
    vector<IJ>::const_iterator itPIJ;
    vector<IJ>::const_iterator itSIJ;
	for (itPIJ = t.Ls.ijs.begin(); itPIJ != t.Ls.ijs.end(); ++itPIJ){
		for (itSIJ = pQ->second.ijs.begin(); itSIJ != pQ->second.ijs.end(); ++itSIJ){
			match_ij(itPIJ, itSIJ, delta_j, pol);
        }
	}
}

Pattern left_convolve(Ptr_pattern_pair pQ, const Pattern& t, const string& prefix_w)
{
	std::vector<IJ> ijs;
	Offset_list ol(ijs);
	Pattern pat(ol);
	Offset_list *pol = &pat.Ls;
    
    //Create new pattern
    string new_pattern = pQ->first + t.motif.substr(prefix_w.size(), string::npos);
    pat.motif = new_pattern;
    //Calculate the offset
    unsigned delta_j = non_bracket_diff_length(pQ->first, prefix_w);
    //Make offset list (pol is pointer to the new pattern's offset list)
    match_ijs_lhs(pQ, t, delta_j, pol);
    
    return pat;
}

Pattern right_convolve(const Pattern& t, Ptr_pattern_pair pQ, const string& suffix_w)
{
	std::vector<IJ> ijs;
	Offset_list ol(ijs);
	Pattern pat(ol);
	Offset_list *pol = &pat.Ls;
	string new_pattern = t.motif + pQ->first.substr(suffix_w.size(), string::npos);
    pat.motif = new_pattern;
    unsigned delta_j = non_bracket_diff_length(t.motif, suffix_w);
    match_ijs_rhs(t, pQ, delta_j, pol);
    
    return pat;
}


// computes the global distance to the IJ
inline unsigned get_global_distance(const IJ& ij, const vector<Sequence>& seqs)
{
	unsigned sum = 0;
    size_t k;
	for (k = 0; k < ij.i; ++k)
		sum += seqs.at(k).str.size();
	return (ij.j + sum);
}

// returns the sum of the distances between each Offset_list element
inline int diff_sum(const vector<IJ>& ijs, const vector<Sequence>& seqs)
{
	int sum = 0;
	size_t len = ijs.size() - 1;
    size_t k;
	for (k = 0; k < len; ++k)
		sum += (get_global_distance(ijs.at(k + 1), seqs) - get_global_distance(ijs.at(k), seqs));
	return sum;
}

Maximal_key make_max_key(const Pattern& pattern, const vector<Sequence>& seqs)
{
	Maximal_key max_key;
    IJ ijpair;
	max_key.offset_list_size = pattern.Ls.ijs.size();
	max_key.diff_sum = diff_sum(pattern.Ls.ijs, seqs);
    max_key.first_diff_sum_sum = get_global_distance(pattern.Ls.ijs.at(1), seqs) - get_global_distance(pattern.Ls.ijs.at(0), seqs);
    max_key.last_diff_sum_sum = get_global_distance(pattern.Ls.ijs.at(pattern.Ls.ijs.size() -1), seqs) - get_global_distance(pattern.Ls.ijs.at(pattern.Ls.ijs.size() - 2), seqs);
    ijpair = pattern.Ls.ijs.front();
    max_key.first_sequence = ijpair.i;
    ijpair = pattern.Ls.ijs.back();
    max_key.last_sequence = ijpair.i;
    
	return max_key;
}


/* Input: Two patterns (string and offset list)                                             */
/* Output: Bool. True, if the second pattern is an offset-wise subpattern of the first.     */
/* False, otherwise.                                                                        */
/* Details: A pattern p2 is not an offset-wise subpattern of p1 iff, they appear on         */
/* different sequences, p1 starts after p2 on the same sequence, p1 ends before p2 on the   */
/* same sequence, and the points from where they start per sequence don't have the exact    */
/* same offset for all appearances                                                          */
bool string_search(const Pattern &max_motif, const Pattern &pat_motif){
    vector<IJ>::const_iterator it1=max_motif.Ls.ijs.begin(), it2=pat_motif.Ls.ijs.begin();
    long int l1 = max_motif.motif.size() - count_bracketed_length(max_motif.motif);
    long int l2 = pat_motif.motif.size() - count_bracketed_length(pat_motif.motif);
    long int offset = it1->j - it2->j;
    //If new is longer literal wise then return false
    while(it1 != max_motif.Ls.ijs.end()){
        if(it1->i != it2->i || it1->j > it2->j || it1->j + l1 < it2->j + l2 || it1->j - it2->j != offset)
            return false;
        it1++;
        it2++;
    }
    return true;
}

//Prints the final results <----------------Check if -v is needed. Technically everything in
//the list satisfies supports requirements.
void out_soln(vector<Pattern>& soln_vec, const Config &config, const vector<Sequence>& seqs)
{
	Soln_Op soln_op;
    unsigned temp;
    vector<IJ>::const_iterator itof;
	sort(soln_vec.begin(), soln_vec.end(), soln_op);
	ofstream cofs(config.output_file.c_str());
	if (cofs.is_open()){
        cofs << "##########################################################" << endl;
        cofs << "#                                                        #" << endl;
        cofs << "#                       FINAL RESULTS                    #" << endl;
        cofs << "#                                                        #" << endl;
        cofs << "##########################################################" << endl;
        vector<Pattern>::const_iterator it;
        vector<IJ>::const_iterator itOL;
        for (it = soln_vec.begin(); it != soln_vec.end(); ++it){
            std::list<unsigned> templist;
            for (itOL = it->Ls.ijs.begin(); itOL != it->Ls.ijs.end(); ++itOL){
                templist.push_back(itOL->i);
            }
            //Finds number per sequence
            templist.unique();
            unsigned temp = templist.size();
            if(!config.k_switch){
                cofs << it->Ls.ijs.size() << "\t" << temp << "\t" << it->motif;
                if (config.offset_report){
                    for(itof = it->Ls.ijs.begin(); itof != it->Ls.ijs.end(); itof++)
                        cofs << " " << itof->i << " " << itof->j;
                }
                cofs << endl;
            }
            else if (temp >= config.k){
                cofs << it->Ls.ijs.size() << "\t" << temp << "\t" << it->motif;
                if (config.offset_report)
                    for(itof = it->Ls.ijs.begin(); itof != it->Ls.ijs.end(); itof++)
                        cofs << " " << itof->i << " " << itof->j;
                cofs << endl;
            }
        }
    }
	else cerr << "Can't open output file." << endl;
}

//Calculates the difference of the actual literal length between
//two strings (a bracket has length 1). To get the literal length
//of just one string, pass the empty string as the second.
unsigned non_bracket_diff_length(string first, string second){
    unsigned length1=0, length2 =0;
    int i=0;
    while (i < first.size()){
        length1++;
        if(first[i] == '[')
            while(first[i] != ']')
                i++;
        i++;
    }
    i=0;
    while (i < second.size()){
        length2++;
        if(second[i] == '[')
            while(second[i] != ']')
                i++;
        i++;
    }
    return length1 - length2;
}

//Returns the kth character of a string.
//Needed when the string includes brackets.
//Starts counting from 0.
int kth_character(int k,string s){
    int i=0,j=0;
    if (k < 0)
        return 0;
    while(j < s.size() && i!=k){
        if(s[j] == '[')
            while (s[j]!=']')
                j++;
        j++;
        i++;
    }
    return j;
}

//Returns the number of brackets in a string
unsigned count_brackets(string s){
    unsigned num = 0;
    int i;
    for(i = 0; i<s.size(); i++)
        if(s[i] == '[')
            num++;
    return num;
}

//Counts the number of characters in brackets (and the brackets) minus 1 per bracket
//length of string - returned value = number of literals
unsigned count_bracketed_length(string s){
    unsigned num =0,i=0;
    while(i<s.size()){
        if(s[i] == '['){
            while(s[i] != ']'){
                num++;
                i++;
            }
        }
        i++;
    }
    return num;
}

/* Input: A string of characters and a character                                            */
/* Output: The string including the new character in alphabetical order                     */
std::string insert_n_order(std::string so, char c){
    int i = 0;
    while(i < so.size() && so[i] < c)
        i++;
    if (i == so.size()){
        so.push_back(c);
        return so;
    }
    string sn = so.substr(0,i);
    sn.push_back(c);
    return sn + so.substr(i);
}


/* Input: The list of offsets of a pattern, the position of the pattern the bracket starts  */
/* and the input file                                                                       */
/* Output: The cleaned up bracket. No non-appearing chars on the bracks body.               */
std::string clean_up_bracket(std::vector<IJ> offsets, long int motif_pos, vector<Sequence> seqs){
    string chars;
    char c;
    vector<IJ>::iterator it;
    for(it = offsets.begin(); it < offsets.end(); it++){
        c = seqs[it->i].str[it->j + motif_pos];
        if (chars.find(c) == std::string::npos)
            chars = insert_n_order(chars,c);
    }
    return chars;
}

/* Input: A pattern and the input file                                                      */
/* Output: The cleaned up motif of a pattern. No brackets have non-appearing characters     */
std::string clean_up_motif(Pattern p, vector<Sequence> seqs){
    long int i=0,m,k=0;
    string s,appeared;
    while(i < p.motif.size()){
        if(p.motif[i] == '['){
            s.push_back('[');
            appeared = clean_up_bracket(p.Ls.ijs, k,seqs);
            for(m=0; m < appeared.size(); m++)
                s.push_back(appeared[m]);
            while (p.motif[i] != ']') i++;
        }
        s.push_back(p.motif[i]);
        i++;
        k++;
    }
    return s;
}

/* Input: All maximal patterns and the input file                                           */
/* Output: All maximal patterns with cleaned up brackets. A clean bracket is one with no    */
/* characters inside that never appeared in the input or for that specific pattern          */
std::vector<Pattern> clean_up_soln(const Maximal_map &maximal, vector<Sequence> seqs)
{
	std::vector<Pattern> solution;
    string s;
    Maximal_map::const_iterator iter;
	for (iter = maximal.begin(); iter != maximal.end(); ++iter){
        s = clean_up_motif(iter->second,seqs);//<----------------ONLY CLEAN BRACKS
        if(s != ""){
            Pattern p(iter->second.Ls.ijs);
            p.motif = s;
            solution.push_back(p);
        }
        //solution.push_back(iter->second); //<--------------ALL BRACKS
    }
	return solution;
}

//Filters the bracket number per sequence
Maximal_map remove_max_bracks(const Maximal_map maximal, unsigned n, vector<Sequence> seqs){
    Maximal_map new_map, new_brack_map;
    Maximal_map::const_iterator iter, iter2;
    unsigned brack_num,i,j;
    vector<int> brack_pos, brack_ends;
    vector<IJ> ijs;
    vector<IJ>::const_iterator iter_ij;
    int pos, pos2, end_pos, start_pos, start_offset;
    string cut_off;
	for (iter = maximal.begin(); iter != maximal.end(); ++iter){
        brack_num = count_brackets(iter->second.motif);
        if (brack_num <= n)
            new_map.insert(std::make_pair(iter->first, iter->second));
        else{
            pos=-1;
            brack_ends.push_back(-1);
            while ((pos = iter->second.motif.find_first_of("[",pos+1)) != string::npos){
                brack_pos.push_back(pos);
                brack_ends.push_back(iter->second.motif.find_first_of("]",pos));
            }
            brack_pos.push_back(iter->second.motif.size());
            for (i=0;i<brack_num-n+1;i++){
                end_pos = brack_pos[i+n]-1;
                start_pos = brack_ends[i] + 1;
                //Find where to start from
                while(iter->second.motif[start_pos] == '.')
                    start_pos++;
                //Count literals in between original starting point and new
                cut_off = iter->second.motif.substr(0,start_pos);
                start_offset = cut_off.size() - count_bracketed_length(cut_off);
                //Change ijs
                for(iter_ij = iter->second.Ls.ijs.begin(); iter_ij != iter->second.Ls.ijs.end(); iter_ij++){
                    IJ ij = {iter_ij->i, iter_ij->j + start_offset};
                    ijs.push_back(ij);
                }
                Offset_list Ls(ijs);
                Pattern p(Ls);
                //Find where to end to
                if(iter->second.motif[end_pos] != iter->second.Ls.ijs.size()-1)
                    while(iter->second.motif[end_pos] == '.')
                        end_pos--;
                p.motif = iter->second.motif.substr(start_pos,end_pos-start_pos+1);
                new_brack_map.insert(std::make_pair(iter->first,p));
                ijs.clear();
            }
            brack_pos.clear();
            brack_ends.clear();
        }
    }
    //If the exact same motif (only motif) exists
    for(iter=new_brack_map.begin(); iter != new_brack_map.end(); ++iter){
        iter2 = new_map.begin();
        while (iter2 != new_map.end() && iter2->second.motif.compare(iter->second.motif)!= 0)
            iter2++;
        if(iter2 == new_map.end())
            new_map.insert(std::make_pair(iter->first,iter->second));
    }
    return new_map;
}

//Filters maximum support per sequence
Maximal_map remove_big_k_patterns(const Maximal_map maximal, unsigned q){
    Maximal_map new_map;
    Maximal_map::const_iterator iter;
    
    for (iter = maximal.begin(); iter != maximal.end(); ++iter){
        if (iter->second.Ls.ijs.size() <= q)
            new_map.insert(std::make_pair(iter->first,iter->second));
    }
    return new_map;
}

