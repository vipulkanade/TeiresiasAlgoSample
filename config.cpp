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

#include "config.h"
#include "pattern.h"

#include <cstring>
#include <iostream>
#include <stdlib.h>

using namespace std;

void print_usage()
{
	cout << endl;
	cout << "======================================================================" << endl;
	cout << "\t\t\tCommand Arguments" << endl;
	cout << "======================================================================" << endl;
	cout << "-h: Turn off help (Command Arguments)." << endl;
	cout << "-w<template_length>, template_length must be unsigned int. \n\t(w < smallest string)" << endl;
	cout << "-l<required_matches>, required_matches must be unsigned_int. (l < w)" << endl;
	cout << "-k<support>, support must be unsigned_int. (k >= 2)" << endl;
	cout << "-i<file_path>, input file_path" << endl;
	cout << "-o<file_path>, default: output.txt" << endl;
	cout << "-v: Requires k support on a per sequence basis." << endl;
	cout << "-b<file_path>, file containing homology list." << endl;
	cout << "-n<unsigned_int>, max_brackets" << endl;
	//cout << "-d: Prints diagnostics." << endl;
	cout << "-p: Prints offset lists." /*<< "\t\t-u: Ignores lowercase."*/ << endl;
	cout << "-s: Scanning phase only." /*<< "\t\t-r: Remove overlaps."*/ << endl;
	cout << "======================================================================" << endl << endl;
}

Config::Config():
	h(true), run(false),
	l(1), w(1), k(1),
	input_file(""),
	output_file(""),
	k_switch(false),
	equivalence_file(""),
	max_brackets(-1),
	/*diagnostics(false),*/
	offset_report(false), /*uppercase(false),*/ /*remove_overlaps(false),*/
    scan_only(false),
    convolution_length(1),
    max_support(-1)
{}

ostream& operator<<(ostream& out, const Config& config)
{
	out << "-h" << config.h << ", run: " << config.run << endl;
	out << "-w" << config.w << ", -l" << config.l << ", -k" << config.k << endl;
	out << "-i" << config.input_file << endl;
	out << "-o" << config.output_file << endl;
	out << "-v: " << config.k_switch << ", -b" << config.equivalence_file <<
		", -n: " << config.max_brackets << endl;
	out /*<< "-d: " << config.diagnostics << ", "*/ << "-p: " << config.offset_report /*<< ", -u: " << config.uppercase*/ <<
		", -s: " << config.scan_only /*<< ", -r: " << config.remove_overlaps*/ << ", -c: " << config.l -1 << " (ALWAYS l-1) "
        ", -q: " << config.max_support << endl;
    out << "Note that -1 for -n and -q means that there is no upper bound for the number of brackets and the support respectively" << endl;
	return out;
}

Config initialize_cmdline(const int argc, const char *const argv[])
{
	Config config;
    int i;
	for (i = 1; i < argc; ++i) // Parse argv[]
	{
		if (!strncmp(argv[i], "-h", 2))
			config.h = false;
		if (!strncmp(argv[i], "-run", 4))
			config.run = true;
		const char *W = "-w";
		int WLen = strlen(W);
		if (!strncmp(W, argv[i], WLen))
			config.w = atoi(&argv[i][WLen]);
		const char *L = "-l";
		int LLen = strlen(L);
		if (!strncmp(L, argv[i], LLen)){
			config.l = atoi(&argv[i][LLen]);
            config.convolution_length = config.l-1;
        }
        const char *Q = "-q";
		int QLen = strlen(Q);
		if (!strncmp(Q, argv[i], QLen))
			config.max_support = atoi(&argv[i][QLen]);
		const char *K = "-k";
		int KLen = strlen(K);
		if (!strncmp(K, argv[i], KLen))
			config.k = atoi(&argv[i][KLen]);
		const char *INPUT_FILE_PATH = "-i";
		const size_t INPUT_VAR_NAME_LEN = strlen(INPUT_FILE_PATH);
		if (!strncmp(INPUT_FILE_PATH, argv[i], INPUT_VAR_NAME_LEN))
			config.input_file.append(&argv[i][INPUT_VAR_NAME_LEN]);
		const char *OUTPUT_FILE_PATH = "-o";
		const size_t OUTPUT_VAR_NAME_LEN = strlen(OUTPUT_FILE_PATH);
		if (!strncmp(OUTPUT_FILE_PATH, argv[i], OUTPUT_VAR_NAME_LEN))
			config.output_file.append(&argv[i][OUTPUT_VAR_NAME_LEN]);
		const char *EQ_FILE_PATH = "-b";
		const size_t EQ_VAR_NAME_LEN = strlen(EQ_FILE_PATH);
		if (!strncmp(EQ_FILE_PATH, argv[i], EQ_VAR_NAME_LEN))
			config.equivalence_file.append(&argv[i][EQ_VAR_NAME_LEN]);
		const char *NUM_BRACKETS = "-n";
		const unsigned NUM_BRACKETS_LEN = strlen(NUM_BRACKETS);
		if (!strncmp(NUM_BRACKETS, argv[i], NUM_BRACKETS_LEN))
			config.max_brackets = atoi(&argv[i][NUM_BRACKETS_LEN]);
		if (!strncmp(argv[i], "-v", 2)) // requires k to be supported on per sequence basis.
			config.k_switch = true;
		/*if (!strncmp(argv[i], "-d", 2))
			config.diagnostics = true;*/
		if (!strncmp(argv[i], "-p", 2))
			config.offset_report = true;
		/*if (!strncmp(argv[i], "-u", 2))
			config.uppercase = true;*/
		if (!strncmp(argv[i], "-s", 2))
			config.scan_only = true;
		/*if (!strncmp(argv[i], "-r", 2))
			config.remove_overlaps = true;*/
	}

	if (config.output_file == "") config.output_file = "output.txt";
	if (config.h)
	{
		print_usage();
		cout << config << endl;
	}
	return config;
}
