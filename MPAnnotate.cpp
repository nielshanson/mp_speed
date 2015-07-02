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

    // Determine database and load matching functional hierarchy information if present
    DB_INFO db_info;
    
    // Detect BLAST.parsed.txt files and information 
    if (options.blast_dir.size() > 0 && options.sample_name.size() > 0) {
        if (options.debug) cout << "\t* Options blast_dir and sample_name options found" << endl;

        getBlastFileNames(options.blast_dir, options.sample_name, options, db_info);

        if (options.debug) {
            cout << "\t* Printing files detected by getBlastFileNames():" << endl;
            cout << "\t* Input_blastouts size: " << db_info.input_blastouts.size() << endl;
            for (unsigned int i =0; i < db_info.input_blastouts.size(); i++ ) {
                cout << "\t\t* db_info: " << db_info.db_names[i] << endl;
                cout << "\t\t* input_blastouts: " << db_info.input_blastouts[i] << endl;
                cout << "\t\t* weight_dbs: " << db_info.weight_dbs[i] << endl;
            }
        }
    }
    
    // Extract and create functional and taxonomic hierachy maps
    vector<string> functional_hierarchy_files = getFunctionalHierarchyFiles(options.functional_categories, options);
    
    // Datastructure to pass database functional hierarchy information to threads
    map<string, map<string, string> > dbNamesToHierarchyIdentifierMaps;
    map<string, IDTREE*> dbNamesToHierarchyIdTree;
    vector<string> db_id_list; // idlist for custom parser for idtree

    string full_path = "";
    string db_name = "";
        
    for ( unsigned int i = 0; i < functional_hierarchy_files.size(); i++ ) {
        full_path = options.functional_categories + "/" + functional_hierarchy_files[i];
        db_name = removeEnding(functional_hierarchy_files[i], ".tree.txt");
        for (unsigned int j = 0; j < db_info.db_names.size(); j++) {
            if ( db_info.db_names[j] == db_name ) {
                dbNamesToHierarchyIdentifierMaps[db_name] = makeHierarchyIdentifierMap(full_path);
                if (db_info.idextractors.find(db_name) == db_info.idextractors.end()) {
                    // Specific parser not found, use genertic idtree parser
                    db_id_list.clear();
                    // Get identifers from database
                    for(map<string, string>::iterator itr = dbNamesToHierarchyIdentifierMaps[db_name].begin();
                        itr != dbNamesToHierarchyIdentifierMaps[db_name].end();
                        itr++) {
                        db_id_list.push_back(itr->first);
                    }
                    // Create idTree
                    IDTREE *idtree = new IDTREE;
                    idtree->createTrie(db_id_list);
                    
                    dbNamesToHierarchyIdTree[db_name] = idtree;
                }
            }
        }
    }
    
    if (options.debug) {
        cout << "Functional hierarchies matched to BLASTout.parsed.txt files:" << endl;
        for(map<string, map<string, string> >::iterator itr = dbNamesToHierarchyIdentifierMaps.begin(); 
            itr != dbNamesToHierarchyIdentifierMaps.end(); 
            itr++) {
            cout << itr->first << endl;
        }
    }
    
    // Load Pathway Tools Trie
    vector<string> ptools_list;
    processPToolsRxnsFile(options.ptools_rxn_file, ptools_list);

    PTOOLS_NODE *root = new PTOOLS_NODE();
    PTOOLS_NODE *start;

    // Iterate through the ptools list of annotations
    char buf[10000]; // Temp buffer
    vector<char *> words; // Vector for parsed fields

    // Parse each line of ptools annotations
    for (std::vector<string>::iterator itr = ptools_list.begin() ; itr != ptools_list.end(); ++itr) {
        split(*itr, words, buf,' ');
        start = root;
        for( unsigned int i=0; i < words.size(); i++) {
            string word = string(words[i]);
            // cout << word << endl;
            if(!start->hasChild(word)) {
                // Create child node
                start->insertChild(word);
                if (i == (words.size()-1)) {
                    // Flag node as finished
                    start->flagChildFinished(word);
                }
            }
            start = start->getChildNode(word);
        }
    }

    if (options.debug) {
        cout << "DEBUG: Loaded tree" << endl;
        // printMetaCycTree(root);
    }

    // Sort the gff file by the orf ids
    if (options.debug) cout << "Sort GFF file by ORF order" << endl;
 
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
        thread_data[i].dbNamesToHierarchyIdTree = dbNamesToHierarchyIdTree;
        thread_data[i].root = root;
    }

    // Create the writer's data object
    WRITER_DATA_ANNOT *writer_data = new WRITER_DATA_ANNOT;
    writer_data->thread_data = thread_data;
    writer_data->num_threads = options.num_threads;
    writer_data->db_info = db_info;
    writer_data->options = options;

    if (options.debug) cout << "Begin processing: " << endl;
    parser.initializeBatchReading();

    while(parser.readBatch()) {  // main loop
        parser.distributeInput(thread_data);
        createThreadsAnnotate(options.num_threads, thread_data, writer_data);
    }
    
    if (options.debug) cout << "Done processing batches" << endl;
    
    // Write out hierarchies
    parser.writeFunctionalHierarchyFiles(writer_data, options);
    
    // write ptools admin files
    writePToolsAdminFiles(options.ptools_dir, options.sample_name);

    if (options.debug) { cout << "End of MPAnnotate()" << endl;}

    return 0;

}
