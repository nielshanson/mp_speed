//
// Created by Niels Hanson on 2015-05-12.
//

#ifndef MP_SPEED_MPCREATEPTOOLSINPUT_H
#define MP_SPEED_MPCREATEPTOOLSINPUT_H

//#include <map>
//#include <algorithm>
//#include <assert.h>
//#include <algorithm>
// #include <vector>
#include <stdio.h>
#include <iostream>
#include <stack>

#include "types.h"
#include "utilities.h"
#include "MPCreatePToolsInputOptions.h"


std::ifstream input;
std::ofstream output;

void processPToolsRxnsFile( string ptools_rxn_file, vector<string> &ptools_list);
void processAnnotationsForPTools(LIST *my_list, PTOOLS_NODE *root, string annotation_file, string ptools_dir);
//string processAnnotationForPtools(vector <char *> annotation_words, PTOOLS_NODE *root);
string processAnnotationForPtools(vector <char *> annotation_words, PTOOLS_NODE *root, PTOOLS_NODE *ptools_ptr,
                                  bool complete, vector <string> word_list, vector <string> max_word_list, string annotation);
void writePfEntry(string orf_id, string annotation_product, string ec_number, int &start_base, int length, std::ofstream &output);
bool pushWordForward(string word, PTOOLS_NODE *& ptools_ptr);
void print_dfs(PTOOLS_NODE *node, string);
void printMetaCycTree(PTOOLS_NODE *root);


#endif //MP_SPEED_MPCREATEPTOOLSINPUT_H
