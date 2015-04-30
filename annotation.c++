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

int processParsedBlastout(string db_name, float weight, string blastoutput, MPAnnotateOptions options, map<string, ANNOTATION> annotation_results) {

    if (options.debug) cout << "In processParsedBlastout()" << endl;

    // Inputs
    string filename = blastoutput; // BLAST/LASTout.parsed.txt
    std::ifstream input; // input filestream
    char buf[1000]; // buffer
    vector <char *> fields; // vector for parsed fields
    int count = 0; // line count

    if (options.debug) {
        cout << "Reading ParsedBlastout: " <<  endl;
        cout << " Filename: " << filename << "\n";
    }

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

/*
 * Given a string with the location of the blast_results directory, getBLASTFileNames returns
 * a list of the out.parsed.txt files.
 */
int getBlastFileNames(string blastdir, string sample_name, DB_INFO &db_info) {
    regex_t regex;
    int reti;
    char msgbuf[1000];
    char tempbuf[1000];

    string algorithm = string("last");

    int i =0;
    for(i =0; i < algorithm.size(); i++ )
        tempbuf[i] = std::toupper(algorithm[i]);
    tempbuf[i]='\0';

    algorithm = string(tempbuf);
    string regexp_text = sample_name + "[.](.*)[.]" + algorithm + "out.parsed.txt";
    reti = compile_regex(&regex, regexp_text.c_str());


    vector<string> files;
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(blastdir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << blastdir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);

    string database;
    cout << "Printing files:" << endl;
    for (unsigned int i = 0;i < files.size();i++) {
        cout << files[i].c_str() << endl;
        database = getpattern(&regex, files[i].c_str(), 1) ;
        if( database.size() > 0)  {
            // cout << database << " " << files[i] <<  endl;
            db_info.db_names.push_back(database);
            db_info.input_blastouts.push_back(files[i]);
            db_info.weight_dbs.push_back(1.0);
        }
    }
    return 0;
}