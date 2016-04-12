#ifndef _alphabet_h
#define _alphabet_h

#include <list>
#include <string>
#include <vector>

using namespace std;


struct char_n_pointers
{
    char alpha;
    std::list<string> inlists;
};

struct charlist
{
    std::vector<char_n_pointers> chars;
};


charlist parse_homology(std::string file_name);
void print_hom (charlist mylist);

#endif