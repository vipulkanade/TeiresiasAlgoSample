#include "config.h"
#include "seqs.h"
#include "pattern.h"
#include "convolution.h"
#include "alphabet.h"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <stack>
#include <string>

using namespace std;

#define LITERAL 1

int main(int argc, char *argv[])
{
	Config config = initialize_cmdline(argc, argv);
	vector<Sequence> seqs = parse_seqs(config.input_file);
	if ((config.l <= config.w) && (config.w <= smallest_seq(seqs)) && (config.l >= 2)) config.run = true;
    
	if (config.run)
	{
		cout << "Scanning..." << endl;
        
		//equivalence file handled by alphabet.h
        charlist subst_alpha;
        if(config.equivalence_file != "")
            subst_alpha = parse_homology(config.equivalence_file);

        //Elementary Patterns are found by applying bit masks on the data
        //Int masks consist of 0 1 2 and bit masks from 0 1
        string bits;
        if(config.equivalence_file == "")
            bits = initialize_bit_mask(config.w, config.l);
        else
            bits = initialize_int_mask(config.w, config.l, config.max_brackets);
        
        //This loop handles the creation of elementary patterns. It generates all permutations of bit masks (sorted high to low)
        //and applies them sequentially on all data. After every application of a bit mask the pattern created is stored in a
        //stack. After the application of a bit mask on a file we count whether each pattern appears more than k times. If it
        //does then we add it to the elementary patterns that are returned by the program
        //vector<Elementary_patterns> dynamic_eps;
        int bracks=0;
        if(config.equivalence_file != ""){
            if(config.l > config.max_brackets)
                bracks = config.max_brackets;
            else
                bracks = config.l;
        }
        Elementary_patterns final_eps;
        Elementary_patterns::iterator eps_it;
        size_t i, j, end;
        string current_bit_mask;
        Elementary_patterns eps;
        do
		{
scan_start:
			current_bit_mask = trim0s(bits);
            //For every line i of the file and for every position j apply the mask and insert the pattern created
            //to a stack
			for (i = 0; i < seqs.size(); ++i){
                end = seqs[i].str.size() - current_bit_mask.size() + 1;
                for (j = 0; j < end; ++j) {
                    make_and_insert(subst_alpha,eps,current_bit_mask,seqs,i,j);
				}
			}
            //Check if patterns on temp list are supported and remove non supported ones from temp list
            eps = supported(eps, config);

            //Add supported ones to the list with the elementary patterns
            eps_it = eps.begin();
            while (eps_it != eps.end()){
                final_eps.insert(*eps_it);
                eps_it++;
            }
            eps.clear();
        }while(prev_permutation(bits.begin() + 1, bits.end()));
        
        //If there id an equivalence file then we need the first bit of the bit mask to be 2 or 1 and NEVER 0 (and only 1
        //if there's no equivalence file). By allowing prev_permutations to change bits from the second position on
        //(that's what the above do-while does) we ensure the never 0 condition. To get both 2 and 1 we start with 2 and
        //go back to start with 1.
        if (bracks > 0 && config.equivalence_file!=""){
            bracks--;
            if (bits[0]=='2' && bracks > 0)
                bits = initialize_int_mask(config.w,config.l,bracks);
            else
                if(bits[0] == '2'){
                    bits = "1" + initialize_int_mask(config.w -1, config.l -1, config.max_brackets);
                    if(config.l > config.max_brackets)
                        bracks = config.max_brackets;
                    else
                        bracks = config.l - 1;

                }
                else{
                    bits = "1" + initialize_int_mask(config.w - 1, config.l-1, bracks);
                }
            goto scan_start;
        }
        
        //If user wants only elementary patterns, print and return
        //Else do the convolution
        vector<Elementary_patterns> dynamic_eps2;
        dynamic_eps2.push_back(final_eps);
        final_eps.clear();
        
        if (config.scan_only)
			out_scan(dynamic_eps2, config, seqs);
		else
		{ //Convolution is done by combining elementary patterns together
            
            cout << "Convolving..." << endl;
            
            //Change file format
			
			const size_t OVERLAP_LEN = config.l - LITERAL;
			Dir_p dir_p = make_dir_p(dynamic_eps2, OVERLAP_LEN);
			Dir_s dir_s = make_dir_s(dynamic_eps2, OVERLAP_LEN);
            
			Maximal_map maximal;
			stack<Pattern> pattern_stack;
            vector<Elementary_patterns>::const_iterator iter;
			Elementary_patterns::const_iterator itP;
            string w, smt;
            unsigned t_bracks, ps_bracks;
            Dir_s::const_iterator itDirS, itDirP;
            vector<Ptr_pattern_pair> u_pepvec;
            vector<Ptr_pattern_pair>::const_iterator itU;
            Ptr_pattern_pair pQ;
            for (iter = dynamic_eps2.begin(); iter != dynamic_eps2.end(); ++iter)
			{ // for...for... == while P not empty from paper (reason: dynamic_ep format)
				for (itP = iter->begin(); itP != iter->end(); ++itP)
				{
                    Pattern p = copy_ep(itP);
					pattern_stack.push(p);
                    while (true)
					{
start:
						if (pattern_stack.empty()){
                            remove_entries(p, dir_p, dir_s, OVERLAP_LEN, '.');
                            break;
                        }
						Pattern t = pattern_stack.top();
                        w = prefix(t.motif, OVERLAP_LEN, '.');
                        t_bracks = count_brackets(t.motif);
                        ps_bracks = count_brackets(w);

                        //Find all patterns whose suffix is the current prefix
                        for(itDirS = dir_s.begin(); itDirS!= dir_s.end(); itDirS++)
                            if(itDirS->first == w){
                                break;
                            }
						if (itDirS != dir_s.end())
						{ // if...for... == while u is not empty from paper (reason: STL way)

							u_pepvec = itDirS->second;
							for (itU = u_pepvec.begin(); itU != u_pepvec.end(); ++itU)
							{
								pQ = *itU; // pointer to minimum element of u_pepvec
                                smt = suffix(pQ->first, OVERLAP_LEN, '.');

                                Pattern r = left_convolve(pQ, t, w);
                                //If new/bigger pattern has the same support as old/smaller one then the new should replace the old
                                //except if we are using brackets. If we have brackets it is possible that two patterns have the same
                                //support size and both must be reported
                                if (config.equivalence_file == "" && r.Ls.ijs.size() == t.Ls.ijs.size())
                                    if (!pattern_stack.empty()) pattern_stack.pop();
                                if ((r.Ls.ijs.size() >= config.k) && is_maximal(maximal, r, seqs))
                                {
                                    pattern_stack.push(r);
                                    goto start;
                                }
							}
						}

						w = suffix(t.motif, OVERLAP_LEN, '.');
                        ps_bracks = count_brackets(w);
                        for(itDirP = dir_p.begin(); itDirP!= dir_p.end(); itDirP++)
                            if(itDirP->first == w){
                                break;
                            }
                        if (itDirP != dir_p.end())
						{ // if...for... == while u is not empty from paper (reason: STL way)
							u_pepvec = itDirP->second;
							for (itU = u_pepvec.begin(); itU != u_pepvec.end(); ++itU)
							{
								pQ = *itU; // pointer to minimum element of u_pepvec
                                smt = prefix(pQ->first, OVERLAP_LEN, '.');
                                Pattern r = right_convolve(t, pQ, w);
                                //If new/bigger pattern has the same support as old/smaller one then the new should replace the old
                                //except if we are using brackets. If we have brackets it is possible that two patterns have the same
                                //support size and both must be reported
                                if (config.equivalence_file == "" && r.Ls.ijs.size() == t.Ls.ijs.size())
                                    if (!pattern_stack.empty()) pattern_stack.pop();
                                
                                //Deciding if we'll report new pattern
                                //If new pattern has support more than k and is maximal then it'll be reported
                                if ((r.Ls.ijs.size() >= config.k) && is_maximal(maximal, r, seqs))
                                {
                                    pattern_stack.push(r);
                                    goto start;
                                }
                            }
						}
						if (pattern_stack.empty()){
                            remove_entries(t, dir_p, dir_s, OVERLAP_LEN, '.');
                            break;
                        }
                        add_pattern(maximal, t, seqs);
						pattern_stack.pop();
					}
				}
			}
            if(config.max_brackets != -1)
                maximal = remove_max_bracks(maximal, config.max_brackets,seqs);
            if(config.max_support != -1)
                maximal = remove_big_k_patterns(maximal, config.max_support);
            vector<Pattern> soln_vec = clean_up_soln(maximal, seqs);
            
			out_soln(soln_vec, config, seqs);
		}
	}
	else cerr << "Is 2 <= L <= W <= minimum sequence length?" << endl;

	return 0;
}
