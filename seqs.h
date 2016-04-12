#ifndef SEQS_H_
#define SEQS_H_

#include "config.h"

#include <string>
#include <vector>

struct Sequence
{
	std::string header;
	std::string str;
};

//Be careful seq_len is real length. It starts from 1 NOT 0.
struct Sequence_line
{
    std::string sequences;
    std::vector<long int> seq_len;
};

// Parses file and returns a vector of sequences
std::vector<Sequence> parse_seqs(const std::string& file_name);

Sequence_line transform_seq (Sequence seqs);

// Error check for min sequence length
size_t smallest_seq(const std::vector<Sequence>& seqs);

// FOR DEBUGGING ONLY: Takes in a vector of sequences and sends them to ostream
void out_seqs(std::ostream& out, const std::vector<Sequence>& seqs);

#endif
