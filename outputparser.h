#ifndef _OUTPUTPARSER
#define _OUTPUTPARSER
#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <map>
#include <regex.h>
#include "utilities.h"
#include "types.h"


using namespace std;


class OutputParser {

private:
        std::ifstream input;
        std::string dbname;
        float ln2;
        float lnk;
        float Lambda;

        int i;
        int BATCH_SIZE; 
        regex_t *r ;
        const char * regex_text;
/*
        self.hits_counts = {}
        self.data = {}
        self.refscores = {}
        self.refBitScores = {}
*/
        map<std::string, bool> query_dictionary;
        map<std::string, unsigned int> refBitScores;
        vector<string> inputbuffer;


public:
      void initialize(); 
      char buf[10000];
      Options options;
      GLOBAL_PARAMS params;
      vector<char *> fields;
      OutputParser(const Options &options, const GLOBAL_PARAMS &params);
      ~OutputParser();
      void create_query_dictionary();
      void create_annotation_dictionary(map<string, string> * annot_map);
      void create_refBitScores(THREAD_DATA *thread_data);
      bool readABatch();
      void closeBatchReading();
      void initializeBatchReading();
      void distributeInput(THREAD_DATA *thread_data);

};


#endif //_OUTPUTPARSER
