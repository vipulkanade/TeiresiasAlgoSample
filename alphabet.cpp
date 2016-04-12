/*****************************************************************************************************/
/* This code was created by Venetia Pliatsika and Jason Mazzatenta and is an implementation of the   */
/* Teiresias algorithm as it appears on the paper Rigoutsos, I, Floratos, A (1998) Combinatorial     */
/* pattern discovery in biological sequences: The TEIRESIAS algorithm. Bioinformatics 14: 55-67.     */
/* Contact us at: https://cm.jefferson.edu/contact-us/                                               */
/*                                                                                                   */
/* Use of these codes is bound by the following terms and conditions:                                */
/*                                                                                                   */
/* Terms of Use:  This code can be freely used for research, academic and other non-profit activities*/
/* (the “Authorized Use”). Commercial use is strictly prohibited.  The code can be copied and        */
/* compiled on any platform for the Authorized Use, but cannot be modified without the written       */
/* permission of the Computational Medicine Center of Thomas Jefferson University                    */
/* https://cm.jefferson.edu                                                                          */
/*                                                                                                   */
/* THE CODE IS PROVIDED “AS IS” WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESSED  */
/* OR IMPLIED. TO THE FULLEST EXTENT PERMISSIBLE PURSUANT TO APPLICABLE LAW. THOMAS JEFFERSON        */
/* UNIVERSITY, AND ITS AFFILIATES, DISCLAIM ALL WARRANTIES, EXPRESS OR IMPLIED, INCLUDING, BUT NOT   */
/* LIMITED TO, THE IMPLIED WARRANTIES OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND*/
/* NON-INFRINGEMENT.                                                                                 */
/*                                                                                                   */
/* NEITHER THOMAS JEFFERSON UNIVERSITY NOR ITS AFFILIATES MAKE ANY REPRESENTATION AS TO THE RESULTS  */
/* TO BE OBTAINED FROM USE OF THE CODE.                                                              */
/*****************************************************************************************************/

#include "alphabet.h"

#include <fstream>
#include <iostream>
#include <list>


using namespace std;

/* Input: The file that describes the substitution rules                                             */
/* A structure that contains pairs of chars that can be substituted and their substitution lists     */
charlist parse_homology(std::string file_name){
    ifstream homf(file_name.c_str());
    charlist mylist;

	if (homf.is_open())
	{
        string line;
        bool flag = false;
        std::vector<char_n_pointers>::iterator it;
        while (homf.good())
		{
            getline(homf, line);
            int i;
            for (i=0; i < line.size(); i++){
                //find if alpha exists in mylist if not add it
                for(it = mylist.chars.begin(); it != mylist.chars.end(); it++){
                    if(it->alpha == line[i]){
                        flag = true;
                        break;
                    }
                }
                //alpha exists
                if(flag){
                    //add reference to the last added element on the list of substitution rules
                    //for this character
                    it->inlists.push_back(line);
                }
                else{ //alpha does not exist / Create new entry
                    char_n_pointers newcharlist;
                    newcharlist.alpha = line[i];
                    newcharlist.inlists.push_back(line);
                    mylist.chars.push_back(newcharlist);
                }
                flag = false;
            }
		}
	} else cerr << "Can't open homology file." << endl;
    
	return mylist;
}

/* Input: The list of characters that can be substituted and the allowed substitutions per character   */
/* Description: Prints on standard output the chars and the substs list for each                       */
void print_hom (charlist mylist){
    std::vector<char_n_pointers>::iterator it;
    for (it=mylist.chars.begin(); it != mylist.chars.end(); ++it){
        cout << it->alpha << " can be replaced by ";
        std::list<string>::iterator it2;
        for(it2=it->inlists.begin(); it2!= it->inlists.end(); it2++){
            cout << *it2 << " ";
        }
        cout << endl;
    }
}