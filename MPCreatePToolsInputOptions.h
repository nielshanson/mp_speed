//
// Created by Niels Hanson on 2015-05-12.
//

#ifndef MP_SPEED_MPCREATEPTOOLSINPUTOPTIONS_H
#define MP_SPEED_MPCREATEPTOOLSINPUTOPTIONS_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <regex.h>

#include "utilities.h"

using namespace std;

// Structure for RPKM input options
struct MPCreatePToolsInputOptions {

    // Processing options
    string ptools_rxn_file;
    string annotation_table;
    string ptools_dir;
    string sample_name;
    bool debug;

    // Constructor with default settings
    MPCreatePToolsInputOptions() {
        ptools_rxn_file = "";
        annotation_table = "";
        ptools_dir = "";
        sample_name = "";
        debug = false;
    }

    bool SetOptions(int argc, char *argv[]);
    void printUsage(char *arg);
    void printOptions();
    bool checkOptions();
    bool set_sample_name_from_annotation_table();
};

#endif //MP_SPEED_MPCREATEPTOOLSINPUTOPTIONS_H
