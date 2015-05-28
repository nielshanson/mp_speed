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
void processPToolsRxnsFile( string ptools_rxn_file, vector<string> &ptools_list);
void processAnnotationsForPTools(LIST *my_list, PTOOLS_NODE *root, string annotation_file);
void processAnnotationForPtools(vector <char *> annotation_words, PTOOLS_NODE *root, LIST *my_list);
bool pushWordForward(char *word, PTOOLS_NODE &ptools_ptr);
void print_dfs(PTOOLS_NODE *node, string);


#endif //MP_SPEED_MPCREATEPTOOLSINPUT_H
