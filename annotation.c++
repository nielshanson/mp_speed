//
// Created by Niels Hanson on 2015-04-29.
//

/*
 * annotation.c++: This class library contains functions associated the processing of multiple parsed database results.
 * (MPAnnotate)
 */

#include "annotation.h"

using namespace std;

#define PRINT_INTERVAL 10000

void readContigLengths(string file, map<string, unsigned int> &contig_lengths) {


    string filename = file;
    std::ifstream input;
    char buf[1000];
    vector <char *> fields;
    int count = 0;

    cout << "Reading refscores " <<  endl;
    cout << "Filename " << filename << "\n";
    input.open(filename.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
        return ;
    }

    // int *counts  = (int *)calloc(options.num_threads, sizeof(int));
    string line;
    string contig_id;

    while( std::getline(input, line ).good()) {
        split(line, fields, buf,'\t');
        if( fields.size()!=3) {
            contig_id.clear();
            input.close();
            return;
        }
        contig_id = fields[0];
        contig_lengths[contig_id] = atoi(fields[2]);



        if (count%PRINT_INTERVAL==0)
            std::cout << "x " << count << std::endl;
        count++;
    }
    input.close();



    std::cout << "Number of contig lengths loaded " <<  count << std::endl;


}

void getBlastFileNames(const MPAnnotateOptions &options, DB_INFO &db_info) {
    // Get blastout files from BlastDB directory

    //TODO: Kishori to write
    // Load dummy annotation results results
    db_info.db_names.push_back("metacyc-v4-2011-07-03");
    db_info.db_names.push_back("COG_2013-12-27");
    db_info.db_names.push_back("CAZY_2014_09_04");

    db_info.input_blastouts.push_back("data/mp_output/hmp_airways_SRS014682/blast_results/hmp_airways_SRS014682.metacyc-v4-2011-07-03.LASTout.parsed.txt");
    db_info.input_blastouts.push_back("data/mp_output/hmp_airways_SRS014682/blast_results/hmp_airways_SRS014682.COG_2013-12-27.LASTout.parsed.txt");
    db_info.input_blastouts.push_back("data/mp_output/hmp_airways_SRS014682/blast_results/hmp_airways_SRS014682.CAZY_2014_09_04.LASTout.parsed.txt");

    db_info.weight_dbs.push_back(1);
    db_info.weight_dbs.push_back(1);
    db_info.weight_dbs.push_back(1);

    return;
}

int processParsedBlastout(string db_name, float weight, string blastoutput, MPAnnotateOptions options, map<string, ANNOTATION> annotation_results) {
    // 'q_length', 'bitscore', 'bsr', 'expect', 'aln_length', 'identity', 'ec'
    // blastparser =  BlastOutputTsvParser(db_name, blastoutput);


    // Inputs
    string filename = blastoutput; // BLAST/LASTout.parsed.txt
    std::ifstream input; // input filestream
    char buf[1000]; // buffer
    vector <char *> fields; // vector for parsed fields
    int count = 0; // line count

    cout << "Reading ParsedBlastout: " <<  endl;
    cout << " Filename: " << filename << "\n";
    input.open(filename.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
        exit(1);
    }

    ANNOTATION annotation;

    string line;
    string query_id; // query_id

    while( std::getline(input, line ).good()) {
        split(line, fields, buf,'\t');
        if( fields.size()!=12) {
            // not a blast file
            query_id.clear();
            input.close();
            continue;
        }
        query_id = fields[0];




        if (count%PRINT_INTERVAL==0)
            std::cout << "x " << count << std::endl;
        count++;
    }
    input.close();



    std::cout << "Number of contig lengths loaded " <<  count << std::endl;
}