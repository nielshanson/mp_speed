//
// Created by Niels Hanson on 2015-04-29.
//

#ifndef _SPEED_ANNOTATION_H
#define _SPEED_ANNOTATION_H

#include <map>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <dirent.h>
#include <stdio.h>
#include <pthread.h>
#include <sstream>

#include "utilities.h"
#include "types.h"
#include "MPAnnotateOptions.h"
#include "idTree.h"

using namespace std;

map<string, string> makeHierarchyIdentifierMap(string hierarchy_filename);

vector<string> getFunctionalHierarchyFiles(string hierarhcy_dir, MPAnnotateOptions options);

HNODE* createHNODE(string line);

// void readContigLengths(string contig_map_file, map<string, unsigned int> &contig_lengths);

int getBlastFileNames(string blastdir, string sample_name, MPAnnotateOptions options, DB_INFO &db_info);

// int processParsedBlastout(string db_name, float weight, string blastoutput, MPAnnotateOptions options, map<string, ANNOTATION> &annotation_results);

string orf_extractor_from_gff(const string &line);

void createThreadsAnnotate(int num_threads, THREAD_DATA_ANNOT *thread_data, WRITER_DATA_ANNOT *writer_data);

void *annotateOrfsForDBs( void *_data) ;

void *reduceAnnotations( void *_writer_data)  ;

ANNOTATION* createAnnotation(const char * line, const string &dbname, bool taxonomy =true);

ANNOTATION* getBestAnnotation(vector<ANNOTATION *> annotation_list, ANNOTATION* annotation);

int process_all_parsed_blastout(THREAD_DATA_ANNOT *thread_data);

string createFunctionalAndTaxonomicTableLine(ANNOTATION annotation);

// void *writeAnnotatedPreamble(void *_writer_data);

void printMetaCycTree(PTOOLS_NODE *root);

void processPToolsRxnsFile( string ptools_rxn_file, vector<string> &ptools_list );

string processAnnotationForPtools(vector <char *> annotation_words, PTOOLS_NODE *root, PTOOLS_NODE *ptools_ptr,
                                  bool complete, vector <string> word_list, vector <string> max_word_list,
                                  string annotation);

bool pushWordForward(string word, PTOOLS_NODE *& ptools_ptr);

void writePfEntry(string orf_id, string annotation_product, string ec_number, string start_base, string end_base, std::ofstream &output);

void writePToolsResults(WRITER_DATA_ANNOT* writer_data, string ptools_dir, string sample_name);

bool  createFunctionWeights(const string &intpufile, const string &outputfile);
#endif //MP_SPEED_ANNOTATION_H
