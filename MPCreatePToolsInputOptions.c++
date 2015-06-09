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
        else if (strncmp(argv[i], "--anno_table", strlen("--anno_table")) == 0) {
            this->annotation_table = argv[++i];
        }
        else if (strncmp(argv[i], "--ptools_dir", strlen("--ptools_dir")) == 0) {
            this->ptools_dir = argv[++i];
        }
        else if (strncmp(argv[i], "--sample_name", strlen("--sample_name")) == 0){
            this->sample_name = argv[++i];
        }
        else if( strncmp(argv[i], "--debug", strlen("--debug")) == 0 ) {
            this->debug = true;
        }
        else {
            cout << "ERROR: Cannot recognize argument " << argv[i] << std::endl;;
            return false;
        }
    } //for loop for arguments processing
    return true;
};

/*
 * Print usage information.
 */
void MPCreatePToolsInputOptions::printUsage(char *arg) {
    std::cout << "USAGE : "   << arg <<"\n"\
              << "      : --ptools_rxns  <metacyc_enzymes_rxns_ecs.txt>                              [REQUIRED]\n"\
              << "      : --anno_table   <functional_and_taxonomic_table.txt>                        [REQUIRED]\n"\
              << "      : --ptools_dir   <ptools_input_dir>                                          [REQUIRED]\n"\
              << "      : --sample_name  <sample_name>                                               [OPTIONAL]\n"\
              << "      : --debug                                                                    [OPTIONAL]\n"
              << std::endl;
};

/*
 * Prints current option settings.
 */
void MPCreatePToolsInputOptions::printOptions() {
        cout << "\nCURRENT OPTIONS:\n"\
         << "---------------\n"
         << "ptools_rxns:    " <<  this->ptools_rxn_file << "\n"\
         << "anno_table:     " <<  this->annotation_table << "\n"\
         << "ptools_dir:     " <<  this->ptools_dir << "\n"\
         << "sample_name:    " <<  this->sample_name << "\n"\
         << "debug:          " <<  this->debug << "\n"
         << "---------------\n\n";
}

/*
 * Checks options while reconfiguring if nessisary
 */

bool MPCreatePToolsInputOptions::checkOptions() {
    if (this->ptools_rxn_file == "") {
        cout << "ERROR: Need to set --ptools_rxns argument to path of metacyc_enzymes_rxns_ecs.txt" << endl;
        return false;
    }
    if (this->annotation_table == "") {
        cout << "ERROR: Need to set --anno_table argument to path of functional_and_taxonomic_table.txt" << endl;
        return false;
    }
    if (this->ptools_dir == "") {
        cout << "ERROR: Need to set --ptools_dir argument to path of ptools/ input directory" << endl;
        return false;
    }
    if (this->sample_name == "") {
        if (!set_sample_name_from_annotation_table()) {
            cout << "ERROR: Could not extract sample_name from functional_and_taxonomic_table.txt, set manually or check other input" << endl;
            return false;
        } else {
            cout << "WARNING: Set sample_name from functional_and_taxonomic_table.txt: " << this->sample_name << endl;
        }
    }
    return true;
}

/*
 * Tries to extract sample_name from set annotation table (i.e., functional_and_taxonomic_table.txt file) which should have
 * form /Users/nielshanson/Dropbox/projects/mp_speed/data/<sample_name>.functional_and_taxonomic_table.txt
 */
bool MPCreatePToolsInputOptions::set_sample_name_from_annotation_table() {
    string temp = this->annotation_table;
    char buf[10000]; // Temp buffer
    std::vector<char *> words; // Vector for parsed fields
    std::vector<char *> words2; // Vector for parsed fields
    split(temp, words, buf,'.');

    // <sample_name> should be first word split by '.'
    temp =  words[0];
    split(temp, words2, buf,'/');
    this->sample_name = words2[words2.size()-1];
    if (this->sample_name != "") {
        return true;
    } else {
        return false;
    }
}