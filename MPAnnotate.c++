//
// Created by Niels Hanson on 2015-04-29.
//

#include "MPAnnotate.h"
#define CHUNK_SIZE 100000

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
    // map<string, unsigned int> contig_lengths; Doesn't seem to be used anymore

    // TODO: contig_lengths does not seem to be used anymore.
    // readContigLengths(options.contig_map_file, contig_lengths);

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

    // Sort the gff file by the orf ids
    if (options.debug) {
        cout << "Sort GFF file" << endl;
    }
    string temp_gff = options.input_gff + ".tmp";
    disk_sort_file(string("/tmp/"), options.input_gff, temp_gff, 1000000, orf_extractor_from_gff);
    remove(options.input_gff.c_str());
    rename(temp_gff.c_str(), options.input_gff.c_str());

    // Initialize MPAnnotateParser
    MPAnnotateParser parser(options, db_info);
    // Create array for thread data
    THREAD_DATA_ANNOT *thread_data = new THREAD_DATA_ANNOT[options.num_threads];

    // Set options for each THREAD_DATA_ANNOT object
    for(unsigned int i = 0; i < options.num_threads; i++) {
        thread_data[i].options = options;
        thread_data[i].db_info = db_info;
    }

    // Create the writer's data object
    WRITER_DATA_ANNOT *writer_data = new WRITER_DATA_ANNOT;

    // Set writer output to array of output streams for each database
    writer_data->output = new std::ofstream[db_info.db_names.size()];

    // Open file streams to temporary annotation_table files
    for(unsigned int i=0; i < db_info.db_names.size(); i++ ) {
        writer_data->output[i].open( db_info.TempFile1(options.output_comp_annot,  db_info.db_names[i]).c_str(), std::ofstream::binary);
    }

    // Update writer
    writer_data->thread_data = thread_data;
    writer_data->num_threads = options.num_threads;
    writer_data->db_info = db_info;

    unsigned int b = 0;
    std::cout << "Begin processing:  \n";
    parser.initializeBatchReading();

    // writeAnnotatedPreamble(writer_data);

    while(parser.readABatch()) {  // main loop
        parser.distributeInput(thread_data);
/*
       for(unsigned int i = 0; i < options.num_threads; i++) 
         thread_data[i].b = b;
*/
        create_threads_annotate(options.num_threads, thread_data, writer_data);
        b = (b+1)%2;
    }
    std::cout << "Done processing batches\n";

    // Sorting the file
    for(unsigned int i =0; i < db_info.db_names.size(); i++ ) {
        std::cout << "Sorting " << db_info.db_names[i] << std::endl;
        disk_sort_file(string("/tmp"),
                       db_info.TempFile1(options.output_gff,  db_info.db_names[i]),
                       db_info.TempFile2(options.output_gff,  db_info.db_names[i]),
                       CHUNK_SIZE,
                       function_extractor_from_list
        );

        db_info.RemoveTempFile1(options.output_gff,  db_info.db_names[i]);

        create_function_weights(
                db_info.TempFile2(options.output_gff,  db_info.db_names[i]),
                db_info.FinalFile(options.output_gff,  db_info.db_names[i])
        );

        db_info.RemoveTempFile2(options.output_gff,  db_info.db_names[i]);

    }

    if (options.debug) { cout << "End of MPAnnotate()" << endl;}

    return 0;

}
