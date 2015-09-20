#ifndef _MPAnnotateParser
#define _MPAnnotateParser
#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <map>
#include <regex.h>
#include <algorithm>

#include "utilities.h"
#include "types.h"
#include "externalsort.h"
#include "MPAnnotateOptions.h"
#include "annotation.h"


using namespace std;

class MPAnnotateParser {

private:
        std::ifstream input;
        std::ofstream output;
        int i;
        DB_INFO db_info;
        int BATCH_SIZE; 
        vector<string> inputbuffer;
        map<string, std::ifstream *> parsed_file_streams;
        map<string, string> leftover_lines;
        map<string, vector<string> > dbwise_inputs; // Database name to database input lines
        map<string, bool> gfforfs;

public:
      void initialize(); 
      char buf[10000];
      MPAnnotateOptions options;
      MPAnnotateParser(const MPAnnotateOptions &options, const DB_INFO & db_info);
      ~MPAnnotateParser();
      bool readBatch();
      void closeBatchReading();
      void initializeBatchReading();
      void distributeInput(THREAD_DATA_ANNOT  *thread_data);
      string prepareRefSeqTaxonomy(string ncbi_id, map<string, string> NCBI_ID_to_Common, NCBITree* ncbi_tree);
      void writeFunctionalHierarchyFiles(WRITER_DATA_ANNOT *writer_data, MPAnnotateOptions options);
};


#endif //_MPAnnotateParser
