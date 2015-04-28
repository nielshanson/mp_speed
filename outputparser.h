#ifndef _OUTPUTPARSER
#define _OUTPUTPARSER
#include <string>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <map>
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
/*
        self.hits_counts = {}
        self.data = {}
        self.refscores = {}
        self.refBitScores = {}
*/
        map<std::string, bool> query_dictionary;
        map<std::string, unsigned int> refBitScores;
        map<std::string, std::string> annot_map;


public:
      void initialize(); 
      char buf[10000];
      Options options;
      GLOBAL_PARAMS params;
      vector<char *> fields;
      OutputParser(const Options &options, const GLOBAL_PARAMS &params);
      ~OutputParser();
      void create_query_dictionary();
      void create_annotation_dictionary();
      void create_refBitScores();

};


#endif //_OUTPUTPARSER
