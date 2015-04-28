#ifndef _UTILITIES
#define _UTILITIES
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <string.h>
#include <stdlib.h>
#include <regex.h>

using namespace std;

// Structure for RPKM input options
struct Options {
   
    string input_blastout;
    string database_name;
    string parsed_output;
    string refscore_file;
    string database_map;
    string algorithm;
    float min_score;
    float min_bsr;
    float max_evalue;
    float lambda, k;
    unsigned int num_threads;
    unsigned int min_length;
    float min_identity;
    unsigned int limit;


    vector<string> read_map_files;

    
    /* Flags and settings */
    bool multi_reads; // flag for detecting multiple mapping of reads
    bool show_status; // shows the counter that countes the number of reads processed, and other info 
                       // on screen
    string reads_map_file_format; // aligner type BWA or BLAST, two SAM files or one
    
    // Constructor with default settings 
    Options(){

       input_blastout ="";
       database_name="";
       parsed_output="";
       refscore_file="";
       database_map ="";
       algorithm ="BLAST";
       min_score = 20;
       max_evalue  = 1e-6;
       min_bsr = 0.4;
       min_length = 30;
       min_identity = 30;
       limit = 5;
       lambda = -1;
       num_threads = 1;
       k = -1;
    };
    
    void print_usage( char *arg);

    void print_options();
 
    bool check_arguments();
    bool SetOptions(int argc, char *argv[]);
	// bool SetOptionCommon( const char *flag, const char *value );
	// bool SetOption( const char *flag, const char *value );
	// bool SetOption2D( const char *flag, const char *value );
	// bool SetOptionEST( const char *flag, const char *value );
	// bool SetOptions( int argc, char *argv[], bool twodata=false, bool est=false );

	void Validate();
	// void ComputeTableLimits( int naan, int typical_len, size_t mem_need );

	void Print();
};

void split(const std::string  &strn, std::vector<char *> &v, char *buf, char d='\t');
std::string get_orf_name(std::string & strn, std::vector<char *> &v, char *buf);


bool matchString(const string &str, const string & stringtomatch, bool fromstart=false);

void get_fasta_sequence_info(const std::string &fasta_file_name);

std::string extract_sequence_name(const std::string &name);

string to_string(unsigned long i);

string ShortenORFId(const string &s, regex_t *r);

string ShortenORFId(const string &s);

#define MAX_ERROR_MSG 0x1000
int compile_regex (regex_t * r, const char * regex_text);

int hashIntoBucket(const char *str, unsigned int index);
#endif //_UTILITIES

