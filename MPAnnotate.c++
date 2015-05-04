//
// Created by Niels Hanson on 2015-04-29.
//

#include "MPAnnotate.h"

using namespace std;

int main( int argc, char** argv) {
    // Parse options

    MPAnnotateOptions options;
    if ( ! options.SetOptions(argc, argv) ) {
        options.printUsage(argv[0]);
        exit(1);
    }

    if (options.debug) { cout << "Beginning of MPAnnotate():" << endl;}

    if (options.debug) {
        options.printOptions();
    }

    // Data structures
    ANNOTATION_RESULTS results_dictionary; // stores results
    map<string, float> dbname_weight; // map to store db_weight TODO: not clear if absolutely needed
    map<string, unsigned int> contig_lengths;

    readContigLengths(options.contig_map_file, contig_lengths);

    cout << options.blast_dir.size() << endl;
    cout << options.sample_name.size() << endl;
    DB_INFO db_info;

    if (options.blast_dir.size() > 0 && options.sample_name.size() > 0) {
        if (options.debug) cout << "\t* Options blast_dir and sample_size options found" << endl;

        getBlastFileNames(options.blast_dir, options.sample_name, options, db_info);

        if (options.debug) {
            cout << "\t* Printing files detected by getBlastFileNames():" << endl;
            cout << "\t* Input_blastouts size: " << db_info.input_blastouts.size() << endl;
            for (unsigned int i =0; i < db_info.input_blastouts.size(); i++ ) {
                cout << "\t\t* db_info: " << db_info.db_names[i] << endl;
                cout << "\t\t* input_blastouts: " << db_info.input_blastouts[i] << endl;
                cout << "\t\t* weight_dbs: " << db_info.weight_dbs[i] << endl;
                dbname_weight[db_info.db_names[i]] = db_info.weight_dbs[i];
            }
        }
    }

    unsigned int priority = 6000;

    // TODO: Sort the gff file by the orf ids
    string temp_gff = options.input_gff + ".tmp";
    disk_sort_file(string("/tmp/"), options.input_gff, temp_gff, 1000000, orf_extractor_from_gff);
    remove(options.input_gff.c_str());
    rename(temp_gff.c_str(), options.input_gff.c_str());



    // Process parsed blastout for each dbname, blastoutput, and weight
    for (unsigned int i=0; i < db_info.input_blastouts.size(); i++ ) {
        int count = processParsedBlastout( db_info.db_names[i],
                                           db_info.weight_dbs[i],
                                           db_info.input_blastouts[i],
                                           options,
                                           results_dictionary[db_info.db_names[i]] );
        priority++; // TODO: Ask Kishori what this priority variable is for
    }

    // Count unique ORFs that got annotation in any database
    map<string, bool> count_annotations;
    ANNOTATION_RESULTS::iterator outer_itr;
    map<string, ANNOTATION>::iterator inner_itr;
    for(outer_itr = results_dictionary.begin(); outer_itr != results_dictionary.end(); outer_itr++) {
        map<string, ANNOTATION> seqname_map = outer_itr->second;
        for (inner_itr = seqname_map.begin(); inner_itr != seqname_map.end(); inner_itr++ ) {
            string seqname = inner_itr->first;
            count_annotations[seqname] = true;
        }
    }

    unsigned int count = count_annotations.size();

    cout << priority << "\tTotal Protein Annotations\t" << count << endl;


    // TODO: Create_annotations from the hits
    createAnnotation(dbname_weight, results_dictionary, options, contig_lengths);


    if (options.debug) { cout << "End of MPAnnotate()" << endl;}

}



