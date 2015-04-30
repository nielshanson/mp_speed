//
// Created by Niels Hanson on 2015-04-29.
//

#ifndef MP_SPEED_ANNOTATION_H
#define MP_SPEED_ANNOTATION_H


#include <map>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include "utilities.h"
#include "types.h"
#include "MPAnnotateOptions.h"

using namespace std;

void readContigLengths(string contig_map_file, map<string, unsigned int> &contig_lengths);

void getBlastFileNames(const MPAnnotateOptions &options, DB_INFO &db_info);

int processParsedBlastout(string db_name, float weight, string blastoutput, MPAnnotateOptions options, map<string, ANNOTATION> annotation_results);

#endif //MP_SPEED_ANNOTATION_H
