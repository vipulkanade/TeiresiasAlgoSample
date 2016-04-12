#include "pattern.h"

#include <fstream>
#include <iostream>
#include <string>

using namespace std;

Offset_list::Offset_list(const vector<IJ>& aijs):
	ijs(aijs)
{}

Pattern::Pattern(const Offset_list& ol):
	Ls(ol)
{}


/**************************         NO BRACKETS       ********************************/
//Initializes a bit mask that will be applied on the data so that we can
//acquire the initial patterns. In the case of no equivalences files we
//use this function that doesn't allow brackets. 1 corresponds to actual
//letters (literals) and 0 to '.' (wildcards).
//If we want brackets we use initialize_int_mask
string initialize_bit_mask(const unsigned w, const unsigned l)
{
	size_t j = w - l;
	string bits;
    size_t i;
	for (i = 0; i < l; ++i)
		bits.push_back('1');
	for (i = 0; i < j; ++i)
		bits.push_back('0');
	return bits;
}

/**************************         BRACKETS         *********************************/
//Initializes a bit mask that will be applied on the data so that we can
//acquire the initial patterns. In the case of equivalences files we use
//this function that allows brackets. 2 are letters (literals) that can
//vary, 1 corresponds to actual specific letters (literals) and 0 to '.'
//(wildcards)
//If we don't want brackets we use initialize_bit_mask
string initialize_int_mask(const unsigned w, const unsigned l, const unsigned max_b){
    size_t c0 = w-l;
    size_t c1;
    string bits;
    unsigned min;
    if (l < max_b){
        min = l;
        c1 = 0;
    }
    else{
        min = max_b;
        c1 = l - max_b;
    }
    size_t i;
    for (i = 0; i < min; ++i)
		bits.push_back('2');
    for (i = 0; i < c1; ++i)
		bits.push_back('1');
	for (i = 0; i < c0; ++i)
		bits.push_back('0');
	return bits;
}

//Applies the bit mask taking into account multiple character appearances in the equivalence file
string make_and_insert(charlist substs, Elementary_patterns& eps, const string& current_bit_mask, const vector<Sequence>& seqs,const unsigned i,const unsigned j)
{
    std::string temp_pattern;
    std::vector<string> pat_to_add;
    std::vector<char_n_pointers>::iterator it;
    std::list<string>::iterator substlists;
    bool flag=false;
    unsigned brack_alter,pat_num,a,b;
    char temp;
    std::string temps;
    
    //Making of the pattern
    pat_to_add.push_back(temp_pattern);
    size_t w;
    for (w = 0; w < current_bit_mask.size(); ++w)
	{
        temp = seqs[i].str[j + w];
        if (current_bit_mask[w] == '1'){
            for(a=0;a<pat_to_add.size();a++)
                pat_to_add[a].push_back(temp);
        }
        else if (current_bit_mask[w] == '0'){
            for(a=0;a<pat_to_add.size();a++)
                pat_to_add[a].push_back('.');
        }
        else if (current_bit_mask[w] == '2'){
            //We need to replace the character with its list of replacements
            //we will replace it with all of them and count the actual occurences
            //on a counting table at the end of the scanning phase we will call a
            //different function and substitute the brackets with the "actually-occured"
            //brackets
            flag = false;
            //Find if the character can be replaced with something
            for (it = substs.chars.begin(); it!= substs.chars.end(); it++){
                if(it->alpha == temp){
                    flag = true;
                    break;
                }
            }
            //There's nothing to replace the character with (no brackets).
            if (!flag){
                return "\0";
            }
            else{
                brack_alter = it->inlists.size() - 1;
                pat_num = pat_to_add.size();
                //add brack_alter - 1 copies of the list in the list. Increase the list by pat_num*(brack_alter-1).
                for(a=0;a<brack_alter;a++)
                    for (b=0;b<pat_num;b++)
                        pat_to_add.push_back(pat_to_add[b]);
                b=0;
                for (substlists = it->inlists.begin(); substlists != it->inlists.end(); substlists++){
                    temps = *(substlists);
                    for(a=0;a<pat_num;a++){
                        pat_to_add[b] = pat_to_add[b] + "[" + temps + "]" + "\0";
                        b++;
                    }
                }
            }
        }
	}
    
    //Insert the elements
    Elementary_patterns::iterator iter;
    vector<IJ> ijs;
    for(a=0;a<pat_to_add.size();a++){
        for(iter = eps.begin(); iter != eps.end(); iter++){
            if(iter->first == pat_to_add[a])
                break;
        }
        //If elementary pattern already exists then just add new occurence to the offset list
        //and increase the counter for the actual characters
        if (iter != eps.end())
        {
            IJ ij = {i, j};
            iter->second.ijs.push_back(ij);
        }
        else
        { //Insert elementary pattern as new
            IJ ij = {i, j};
            ijs.push_back(ij);
            Offset_list Ls(ijs);
            eps.insert(make_pair(pat_to_add[a], Ls));
            ijs.clear();
        }
    }
    
    return "\0";
}


//Checking if the patterns in a temporary list are supported, that is if each pattern appears more
//k (input parameter) per input line or in general, in the entire input.
Elementary_patterns supported(Elementary_patterns& eps, const Config& config) {
    Elementary_patterns new_eps;
    //Checking appearances for the entire input
    Elementary_patterns::const_iterator iter;
    for (iter = eps.begin(); iter != eps.end(); ++iter)
        if (iter->second.ijs.size() >= config.k)
            new_eps.insert(*iter);
    return new_eps;
}


inline void out_eps(ofstream& cofs, const Elementary_patterns& eps, const Config& config, const vector<Sequence>& seqs)
{
    Elementary_patterns::const_iterator iter;
	std::list<unsigned> templist;
    vector<IJ>::const_iterator itOL;
    vector<IJ>::const_iterator iterIJ;
    for (iter = eps.begin(); iter != eps.end(); ++iter)
	{
        for (itOL = iter->second.ijs.begin(); itOL != iter->second.ijs.end(); ++itOL){
            templist.push_back(itOL->i);
        }
        templist.unique();
        if (!config.k_switch){
            cofs << iter->second.ijs.size() << "\t" << templist.size() << "\t" << iter->first;
            if (config.offset_report)
            {
                for (iterIJ = iter->second.ijs.begin(); iterIJ != iter->second.ijs.end(); ++ iterIJ)
                    cofs << " " << iterIJ->i << " " << iterIJ->j;
            }
            cofs << endl;
        }
        else if (templist.size() >= config.k){
            cofs << iter->second.ijs.size() << "\t" << templist.size() << "\t" << iter->first;
            if (config.offset_report)
            {
                for (iterIJ = iter->second.ijs.begin(); iterIJ != iter->second.ijs.end(); ++ iterIJ)
                    cofs << " " << iterIJ->i << " " << iterIJ->j;
            }
            cofs << endl;
        }
        templist.clear();
	}
}

void out_scan(const vector<Elementary_patterns>& dynamic_eps, const Config& config, const vector<Sequence>& seqs)
{
	ofstream cofs(config.output_file.c_str());
	if (cofs.is_open()){
        vector<Elementary_patterns>::const_iterator iter;
        for (iter = dynamic_eps.begin(); iter != dynamic_eps.end(); ++iter)
		{
			out_eps(cofs, *iter, config, seqs);
		}
    }
	else cerr << "Can't open output file." << endl;
}


//This function erases all 0s at the end of a string
string trim0s(const string& bit_mask)
{
	string new_bit_mask = bit_mask;
	string::reverse_iterator rit = new_bit_mask.rbegin();
	while(*rit++ == '0') // trim trailing
		new_bit_mask.erase(new_bit_mask.end() - 1);
	return new_bit_mask;
}

