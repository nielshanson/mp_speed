//
// Created by Niels Hanson on 2015-04-29.
//

#ifndef MP_SPEED_MPANNOTATEOPTIONS_H
#define MP_SPEED_MPANNOTATEOPTIONS_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <regex.h>

using namespace std;

// Structure for RPKM input options
struct MPAnnotateOptions {

    // Processing options
    string input_blastout;
    string algorithm;
    string contig_map_file;
    string database_name;
    string blast_dir;
    string rRNA_16S;
    string tRNA;
    bool taxonomy;
    string output_gff;
    string output_comp_annot;
    string input_gff;
    string sample_name;
    float weight_db;

    // Threading options
    unsigned int num_threads;

    // DEBUG options
    bool debug;

    // Constructor with default settings
    MPAnnotateOptions() {
        input_blastout = "";
        algorithm = "BLAST";
        contig_map_file = "";
        database_name = "";
        blast_dir = "";
        rRNA_16S = "";
        tRNA = "";
        sample_name = "";
        taxonomy = false;
        output_gff = "";
        output_comp_annot = "";
        input_gff = "";
        num_threads = 1;
        weight_db = 1.0;
        debug = false;
    };

    void printUsage(char *arg);
    void printOptions();

    // bool CheckArguments();
    bool SetOptions(int argc, char *argv[]);

};

#endif //MP_SPEED_MPANNOTATEOPTIONS_H
