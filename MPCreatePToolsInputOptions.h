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

using namespace std;

// Structure for RPKM input options
struct MPCreatePToolsInputOptions {

    // Processing options
    string ptools_rxn_file;
    // Constructor with default settings
    MPCreatePToolsInputOptions() {
        ptools_rxn_file = "";
    }

    bool SetOptions(int argc, char *argv[]);
};

#endif //MP_SPEED_MPCREATEPTOOLSINPUTOPTIONS_H
