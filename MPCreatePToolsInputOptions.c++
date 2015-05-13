//
// Created by Niels Hanson on 2015-05-12.
//

#include "MPCreatePToolsInputOptions.h"

/*
 * Parse arguments for MetaPathways annotation step.
 */
bool MPCreatePToolsInputOptions::SetOptions(int argc, char *argv[]) {
    for(int i = 1; i < argc ; i++) {
        if (strncmp(argv[i], "--ptools_rxns", strlen("--ptools_rxns")) == 0) {
            this->ptools_rxn_file = argv[++i];
        }
        else {
            cout << "ERROR: Cannot recognize argument " << argv[i] << std::endl;;
            return false;
        }
    } //for loop for arguments processing

    return true;
};