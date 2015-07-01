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
    string algorithm;
    string blast_dir; // blast_results folder
    string functional_categories; // location of functional hierarchy .tree.txt files
    string sample_name;
    string input_gff;
    bool taxonomy;
    string results_dir;
    float weight_db;

    // Threading options
    unsigned int num_threads;

    // DEBUG options
    bool debug;
    
    // Constructor with default settings
    MPAnnotateOptions() {
        algorithm = "BLAST";
        blast_dir = "";
        functional_categories = "";
        sample_name = "";
        input_gff = "";
        taxonomy = true;
        results_dir = "";
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
