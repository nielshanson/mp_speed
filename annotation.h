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

#include "utilities.h"
#include "types.h"
#include "MPAnnotateOptions.h"

using namespace std;

void readContigLengths(string contig_map_file, map<string, unsigned int> &contig_lengths);

int getBlastFileNames(string blastdir, string sample_name, MPAnnotateOptions options, DB_INFO &db_info);

int processParsedBlastout(string db_name, float weight, string blastoutput, MPAnnotateOptions options, map<string, ANNOTATION> &annotation_results);

void createAnnotation(map<string, float> dbname_weight, ANNOTATION_RESULTS results_dictionary, MPAnnotateOptions options, map<string, unsigned int> contig_lengths);

string orf_extractor_from_gff(const string &line);

void createThreadsAnnotate(int num_threads, THREAD_DATA_ANNOT *thread_data, WRITER_DATA_ANNOT *writer_data);

void *annotateOrfsForDBs( void *_data) ;

void *writeAnnotatedGFFs( void *_data) ;

ANNOTATION * createAnnotation(const char * line, const string &dbname, bool taxonomy =true);

int process_all_parsed_blastout(THREAD_DATA_ANNOT *thread_data);

// void *writeAnnotatedPreamble(void *_writer_data);

bool  createFunctionWeights(const string &intpufile, const string &outputfile);
#endif //MP_SPEED_ANNOTATION_H
