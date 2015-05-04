#ifndef __PARSE_TYPE
#define __PARSE_TYPE
#include <utility>
#include <vector>
#include <iostream>
#include <map>

#include "options.h"
#include "MPAnnotateOptions.h"
#include "utilities.h"

typedef struct _GLOBAL_PARAMS {
   float lambda;
   float k;
} GLOBAL_PARAMS;


typedef struct _BLASTOUT_DATA BLASTOUT_DATA;

typedef struct _THREAD_DATA {
   vector<string>  lines;
   vector<BLASTOUT_DATA>  output[2];
   vector< std::pair<std::string, std::string> > refscorePairs;
   map<std::string, unsigned int> refscores;
   float ln2, lnk, lambda;
   Options options;
   map<std::string, std::string> *annot_map;
   unsigned int b;
} THREAD_DATA;


typedef struct _DB_INFO {
    vector<string> db_names;
    vector<string> input_blastouts;
    vector<float> weight_dbs;
} DB_INFO;





typedef struct _BLASTOUT_DATA {
    string query, target, product, ec, taxonomy; 
    unsigned int q_length, aln_length ;
    float bitscore, bsr, expect, identity;
    void setDefault() {
         query = "";
         target = "";
         product = "";
         ec="";
         
         q_length =0;
         aln_length=0;
         bitscore=0;
         bsr=0;
         expect=100;
         identity =0;
    }
} BLASTOUT_DATA;

typedef struct _WRITER_DATA {
    THREAD_DATA *thread_data;
    unsigned int num_threads;
    std::ofstream output;
} WRITER_DATA;


// Annotation data structure
typedef struct _ANNOTATION {
    float bsr, value; 
    string ec, product, taxonomy;
} ANNOTATION;

typedef map<string, map<string, ANNOTATION> > ANNOTATION_RESULTS;

typedef struct _THREAD_DATA_ANNOT {
   vector<ANNOTATION>  lines;

   ANNOTATION_RESULTS annot_objects;

   vector<ANNOTATION>  output[2];
   MPAnnotateOptions options;
   unsigned int b;
   unsigned int tot_annot;
   map<string, unsigned int> dbwise_annot_count;
   DB_INFO db_info; 

   _THREAD_DATA_ANNOT() {
     tot_annot=0;
   }
} THREAD_DATA_ANNOT;

typedef struct _WRITER_DATA_ANNOT {
    THREAD_DATA_ANNOT  *thread_data;
    unsigned int num_threads;
    std::ofstream output;
} WRITER_DATA_ANNOT;

// Annotation map
typedef map<string, map<string, ANNOTATION> > ANNOTATION_RESULTS;

#endif //__RPKM_TYPE



