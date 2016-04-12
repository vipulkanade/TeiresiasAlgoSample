#ifndef CONFIG_H_
#define CONFIG_H_

#include <string>

// TODO: implement k_switch, diagnostics, uppercase, equivalence_file, max_brackets, and remove_overlaps
struct Config
{
	Config();
	bool run; // execute program
	bool h; // get help (command argument list)(Defaults to true, user must turn it off to not display it)
	unsigned l; // at least this many residues (l <= w)
	unsigned w; // elementary pattern (P) size
	unsigned k; // min support (k <= number of seqs in input file)(has two forms and switch, see k_switch)
	std::string input_file;
	std::string output_file;
	bool k_switch; // switch from min k per sequence to min k in general (found anywhere). Default = false means general
	std::string equivalence_file; // homology list. "" means that there are no homologies to be checked
	int max_brackets; // max number of homologies (default = -1)
	//bool diagnostics; // if true, provides diagnostics for output
	bool offset_report; // if true, reports {seq_id, offset} info for patterns
	//bool uppercase; // if true, ignores any patterns that contain lowercase ASCII characters
	bool scan_only; // perform scan only and report Ps
    int convolution_length; //inactive and set up to be l-1 always
    int max_support; //maximum support of pattern
	//remove_overlaps;
};

// Initialize the command line arguments
Config initialize_cmdline(const int argc, const char *const argv[]);

#endif
