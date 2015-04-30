//
// Created by Niels Hanson on 2015-04-29.
//

#include "MPAnnotate.h"

using namespace std;

int main( int argc, char** argv) {
    // Parse options
    cout << "Hello from MPAnnotate" << endl;
    MPAnnotateOptions options;
    if ( ! options.SetOptions(argc, argv) ) {
        options.print_usage(argv[0]);
        exit(1);
    }

    // Data structures
    ANNOTATION_RESULTS results_dictionary; // stores results
    map<string, float> db_name_weight;
    map<string, unsigned int> contig_lengths;

    readContigLengths(options.contig_map_file, contig_lengths);

    cout << options.blast_dir.size() << endl;
    cout << options.sample_name.size() << endl;

    DB_INFO db_info;
    if (options.blast_dir.size() > 0 && options.sample_name.size() > 0) {
        cout << "found" << endl;

        getBlastFileNames(options, db_info);
    } else {
        db_info.db_names.push_back(options.database_name);
        db_info.input_blastouts.push_back(options.sample_name);
        db_info.weight_dbs.push_back(options.weight_db);
    }

    const unsigned int priority = 6000;

    // Process parsed blastout for each dbname, blastoutput, and weight
    for (int i =0; i < db_info.input_blastouts.size(); i++ ){
        int count = processParsedBlastout( db_info.db_names[i], db_info.weight_dbs[i], db_info.input_blastouts[i], options, results_dictionary[db_info.db_names[i]] );
    }




    map<string, bool> count_annotations;
    // TODO: Count unique ORFs that got annotation in any database
    //    for dbname in results_dictionary:
    //        for seqname in results_dictionary[dbname]:
    //            count_annotations[seqname] = True
    //            count = len(count_annotations)
    //    if runstatslogger!=None:
    //        runstatslogger.write("%s\tTotal Protein Annotations\t%s\n" %( str(priority),  str(count)))

    // TODO: Create_annotations from the hits
    // create_annotation(dbname_weight, results_dictionary, opts.input_gff, opts.rRNA_16S, opts.tRNA, opts.output_gff, opts.output_comparative_annotation, contig_lengths)





    cout << "Hello from MPAnnotate" << endl;

}



