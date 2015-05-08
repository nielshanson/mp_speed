#ifndef _MPAnnotateParser
#define _MPAnnotateParser
#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <map>
#include <regex.h>

#include "utilities.h"
#include "types.h"
#include "externalsort.h"
#include "MPAnnotateOptions.h"
#include "annotation.h"


using namespace std;

class MPAnnotateParser {

private:
        std::ifstream input;
       
        int i;
        DB_INFO db_info;
        int BATCH_SIZE; 
        vector<string> inputbuffer;
        map<string, std::ifstream *> parsed_file_streams;
        map<string, string> leftover_lines;
        map<string, vector<string> > dbwise_inputs;
        map<string, bool>  gfforfs;



public:
      void initialize(); 
      char buf[10000];
      MPAnnotateOptions options;
      MPAnnotateParser(const MPAnnotateOptions &options, const DB_INFO & db_info);
      ~MPAnnotateParser();
      bool readABatch();
      void closeBatchReading();
      void initializeBatchReading();
      void distributeInput(THREAD_DATA_ANNOT  *thread_data);
};


#endif //_MPAnnotateParser