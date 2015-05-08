#ifndef __PARSE_TYPE
#define __PARSE_TYPE
#include <utility>
#include <vector>
#include <iostream>
#include <map>
#include <stdio.h>

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


typedef enum {KEGG, COG, SEED } DBTYPE;

typedef struct _DB_INFO {
    vector<string> db_names;
    vector<string> db_tags;
    vector<string> input_blastouts;
    vector<float> weight_dbs;
    vector<DBTYPE> dbtypes;
    vector<string (*)(const char *)> idextractors;

    string TempFile1(const string &str, const string &db)  {
      string name;
      name = str + string(".") + db+ string(".")+ string("DB_annotation_table.txt.tmp1");
      return name;
    }
    
    string TempFile2(const string &str, const string &db)  {
      string name;
      name = str + string(".") + db+ string(".")+ string("DB_annotation_table.txt.tmp2");
      return name;
    }
    
    string FinalFile(const string &str, const string &db)  {
      string name;
      name = str + string(".") + db+ string(".")+ string("DB_annotation_table.txt");
      return name;
    }

    void RemoveTempFile1(const string &str, const string &db)  {
      string name;
      name = str + string(".") + db+ string(".")+ string("DB_annotation_table.txt.tmp1");
      remove(name.c_str());
    }

    void  RemoveTempFile2(const string &str, const string &db)  {
      string name;
      name = str + string(".") + db+ string(".")+ string("DB_annotation_table.txt.tmp2");
      remove(name.c_str());
    }


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

typedef map<string, map<string, ANNOTATION *> > ANNOTATION_RESULTS;

typedef vector<string> DB_HIT;

typedef struct _THREAD_DATA_ANNOT {
   vector<ANNOTATION>  lines;
   ANNOTATION_RESULTS annot_objects;
   vector<ANNOTATION>  output[2];
   vector<string> orfids;
   vector<short int> annot_from_db;
   vector<DB_HIT *> db_hits;

   MPAnnotateOptions options;
   unsigned int b;
   unsigned int tot_annot;
   map<string, unsigned int> dbwise_annot_count;
   DB_INFO db_info; 

   _THREAD_DATA_ANNOT() {
     tot_annot=0;
   }

   void clear() {
      lines.clear();

//      annot_object.clear();
      orfids.clear();
 //     annot_from_db.clear();
  //    db_hits.clear();


   }

} THREAD_DATA_ANNOT;


typedef struct _WRITER_DATA_ANNOT {
    THREAD_DATA_ANNOT  *thread_data;
    unsigned int num_threads;
    std::ofstream *output;
    DB_INFO db_info; 
} WRITER_DATA_ANNOT;

// Annotation map
typedef map<string, map<string, ANNOTATION *> > ANNOTATION_RESULTS;

#endif //__RPKM_TYPE



