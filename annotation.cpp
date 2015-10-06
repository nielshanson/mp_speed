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

//void readContigLengths(string file, map<string, unsigned int> &contig_lengths) {
//
//    string filename = file;
//    std::ifstream input;
//    char buf[1000];
//    vector <char *> fields;
//    int count = 0;
//
//    input.open(filename.c_str(), std::ifstream::in);
//    if(!input.good()){
//        std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
//        return ;
//    }
//
//    // int *counts  = (int *)calloc(options.num_threads, sizeof(int));
//    string line;
//    string contig_id;
//
//    while( std::getline(input, line ).good()) {
//        split(line, fields, buf,'\t');
//        if( fields.size()!=3) {
//            contig_id.clear();
//            input.close();
//            return;
//        }
//
//        contig_id = fields[0];
//        contig_lengths[contig_id] = atoi(fields[2]);
//
//        if (count%PRINT_INTERVAL==0)
//            std::cout << "x " << count << std::endl;
//        count++;
//    }
//    input.close();
//
//    std::cout << "Number of contig lengths loaded " <<  count << std::endl;
//
//}

/*
 * Takes a functional hierarchy file and returns a map<string, string> of hierarchy identifiers to 
 * aliases (if present).
 */
map<string, string> makeHierarchyIdentifierMap(string hierarchy_filename) {
    // cout << "In makeHierarchyIdentifierMap()" << endl;
    
    map<string, string> h_map;  // Hierarchy map
    HNODE *pnode = new HNODE;
    HNODE *cnode = new HNODE;
    // vector<HNODE *> node_stack;
    
    pnode->depth = -1;
    pnode->name = "root";
    pnode->alias = "root";
    
    string line;
    std::ifstream input;
    input.open(hierarchy_filename.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cout << "Error opening '"<< hierarchy_filename <<"'" << std::endl;
        exit(-1);
    }
    
    while( std::getline(input, line ).good()) {
        cnode = createHNODE(line); // get current node
        if (cnode->depth <= pnode->depth) {
            // current node on same level or higher
            h_map[pnode->name] = pnode->alias;
        }
        pnode = cnode;
    }
    
    // add final line
    h_map[pnode->name] = pnode->alias;
    
    input.close();
    
    // // Print out map values
    // for(map<string, string>::iterator itr = h_map.begin(); itr != h_map.end(); ++itr) {
    //     cout << itr->first << "\t" << itr->second << endl;
    // }
    
    return h_map;
}

/* 
 * Given a line from a functional_category *tree.txt file, this function calculates the node's
 * depth in the tree and parses the database id and alias. Returns a pointer to the created HNODE.
 */ 
HNODE* createHNODE(string line) {
    
    HNODE *new_node = new HNODE;
    
    // Count depth
    int depth = 0;
    string name = "";
    string alias = "";
    for (unsigned int i=0; i < line.size(); i++) {
        if (line[i] == '\t') {
            depth++;
        }
    }
    
    // Split line by tabs
    char buf[10000]; // TODO: Check to see if this buffer is needed.
    vector <char *> words;
    split(line, words, buf,'\t');
    
    // parse out name and alias fields
    bool first = true;
    for (unsigned int i=0; i < words.size(); i++) {
        if (hasCharacter(words[i])) {
            if (first) {
                name = words[i];
                first = false;
            } else {
                alias = words[i];
            }
        }
    }
    
    // Fill in new node
    new_node->name = name;
    new_node->alias = alias;
    new_node->depth = depth;
    
    return(new_node);
}

/*
 * Given a product string this function computes the products word information score. I.e., one point for every
 * non-trivial descripive word.
 */
float wordInformation(const string &product) {

    // Split product into individual words
    char buf[10000]; // TODO: Check to see if this buffer is needed.
    vector <char *> words;
    split(product, words, buf,' ');

    float word_information_score = 0.0; // information score of the current words

    // Prepare map of stop words to check for
    static const string arr[] = {"", "is", "have", "has", "will", "can", "should",  "in",
                                 "at", "upon", "the", "a", "an", "on", "for", "of", "by",
                                 "with" ,"and",  ">", "predicted", "protein", "conserved" };

    map <string, int> stop_words;
    for (unsigned int i = 0; i < arr->size(); i++) {
        stop_words[arr[i]] = 1;
    }

    // Check each word for membership in stop_words and presence of underscores '_'
    for (unsigned int i = 0; i < words.size(); i++) {
        string my_word = words[i];
        if (stop_words.count(my_word) <= 0) {
            // not a stop word
            if (my_word.find("_") == std::string::npos) {
                // does not contain an underscore
                word_information_score += 1.0;
            }
        }
    }

    return word_information_score;
}

/*
 * Computes the annotation value of the given ANNOTATION object. Returns an annotation value based on the presence
 * of Enzyme Commission (EC) numbers (+10) and the number of non-trivial words in the ANNOTATION product field.
 */
float computeAnnotationValue(ANNOTATION *annotation) {

    float score = 0.0; // overall annotation score

    if (annotation->ec != "") {
        // annotation object has non-trivial EC number
        score += 10;
    }

    if (annotation->product.find("hypothetical protein") != std::string::npos) {
        score += wordInformation(annotation->product);
    }

    return score;

}

/*
 * Given a B/LASTout.parsed.txt line, this function parses fields and creates an ANNOTATION
 * object returning a pointer to that object.
 */ 
ANNOTATION* createAnnotation(const char * line, const string &dbname, bool taxonomy) {
     vector<char *>fields;
     char buf[10000];  
     string db_name;
     ANNOTATION *annotation = new ANNOTATION;

     split(line, fields, buf,'\t');
     if( fields.size()!=10) {
         return 0;
     }

     try {
        annotation->query = string(fields[0]);
        annotation->target = string(fields[1]);
        annotation->bitscore =  atof(fields[3]);
        annotation->bsr = atof(fields[4]);
        annotation->ec = string(fields[8]);
        annotation->product = string(fields[9]);
        annotation->length = atoi(fields[6]);
        db_name = to_upper(dbname);

        if (taxonomy && (db_name.find("REFSEQ") != std::string::npos)) {
             // TODO: Could be more reliable to get taxonomy from GI number.
             annotation->taxonomy = getTaxonomyFromProduct(annotation->product.c_str());
        }

    }
    catch(...) {
        return 0;
    }

    return annotation;
}

//int processParsedBlastout(string db_name, float weight, string blastoutput, MPAnnotateOptions options, map<string, ANNOTATION> &annotation_results) {
//
//    if (options.debug) cout << "In processParsedBlastout()" << endl;
//
//    // Prepare inputs and buffers
//    string filename = options.blast_dir + "/" + blastoutput; // BLAST/LASTout.parsed.txt
//    std::ifstream input; // Input filestream
//    char buf[10000]; // Temp buffer TODO check to see if multiple buffers nessisary
//    vector <char *> fields; // Vector for parsed fields
//
//    int count = 0; // Line count
//
//    // char tempbuf[1000];
//
//    if (options.debug) {
//        cout << "Reading ParsedBlastout: " << filename << "\n";
//    }
//
//    // Open the file.
//    input.open(filename.c_str(), std::ifstream::in);
//    if(!input.good()){
//        std::cerr << "Error opening '" << filename << "'. Bailing out." << std::endl;
//        exit(1);
//    }
//
//    // ANNOTATION fields for parsing BLAST/LASTout.annotated.txt files
//    ANNOTATION annotation;
//    string line;
//    string query_id;
//    string bsr;
//    string ec;
//    string product;
//    string taxonomy = "";
//
//    // Parse annotation hit line, create ANNOTATION objects
//    while( std::getline(input, line ).good()) {
//        split(line, fields, buf,'\t');
//        if (count == 0) { count++; continue; };
//
//        if( fields.size()!=10) {
//            // Not a parsed blast file
//            cerr << "Parsed BLAST/LASTout file " <<  filename << " did not have the 10 columns." << endl;
//            query_id.clear();
//            input.close();
//            return 1;
//        }
//
//        query_id = fields[0];
//        bsr = fields[4];
//        ec = fields[8];
//        product = fields[9];
//
//        db_name = to_upper(db_name);
//
//        // Construct annotation
//        ANNOTATION annotation;
//        annotation.bsr = atof(bsr.c_str());
//        annotation.ec = ec;
//        annotation.product = product;
//
//        // Extract RefSeq taxonomy from product field
//        if (options.taxonomy && (db_name.find("REFSEQ") != std::string::npos)) {
//            // TODO: Could be more reliable to get taxonomy from GI number.
//            taxonomy = getTaxonomyFromProduct(product.c_str());
//        }
//        annotation.taxonomy = taxonomy;
//
//        // Compute information score of current annotation.
//        annotation.value = computeAnnotationValue(&annotation) * weight;
//
//        // Check to see if current query_id is already in this database's annotation_results
//        if (annotation_results.count(query_id) <= 0) {
//            annotation_results[query_id] = annotation;
//        }
//        else {
//            // Replace existing ANNOTATION with the current annotation if strong annotation score
//            if (annotation_results[query_id].value < annotation.value) {
//                annotation_results[query_id] = annotation;
//            }
//        }
//
//        if (options.debug) {
//            if (count%PRINT_INTERVAL==0)
//                std::cout << "x " << count << std::endl;
//        }
//
//        count++;
//    }
//    input.close();
//
//    if (options.debug) {
//        std::cout << "Number of parsed BLAST/LAST results loaded " <<  count << std::endl;
//    }
//    return count;
//}


/*
 * Given a string with the location of functional_categories directory, getFunctionalHierarchyFiles returns
 * a list of the <database_name>.tree.txt files in the  MetaPathwaysDBs/functional_categories/
 */ 
vector<string> getFunctionalHierarchyFiles(string hierarhcy_dir, MPAnnotateOptions options) {
    
    string hier_ending = ".tree.txt";
    
    // Open directory
    vector<string> files; // list of files in directory
    DIR *dp; // directory pointer
    struct dirent *dirp;
    if((dp  = opendir(hierarhcy_dir.c_str())) == NULL) {
        cout << "Error: opening directory " << hierarhcy_dir << endl;
    }
    
    // Add files in directory to files vector
    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    
    vector<string> hier_files;
    
    // Check to see if files have .tree.txt ending
    for (unsigned int i = 0; i < files.size(); i++) {
        if (hasEnding(files[i], hier_ending)) {
            hier_files.push_back(files[i]);
        }
    }
    
    return(hier_files);
    
}

/*
 * Given a string with the location of the blast_results directory, getBLASTFileNames returns
 * a list of the <DB>.B/LASTout.parsed.txt files.
 */
int getBlastFileNames(string blastdir, string sample_name, MPAnnotateOptions options, DB_INFO &db_info) {

    regex_t regex; // Regular expression.
    // char tempbuf[1000]; // Buffer.

    // Convert algorithm name to uppercase (i.e., BLAST, LAST)

    string algorithm = options.algorithm;

    algorithm = to_upper(algorithm);

    // Create regular expression pattern for BLAST/LASTout.parsed.txt files
    string regexp_text = sample_name + "[.](.*)[.]" + algorithm + "out.parsed.txt$";
    compile_regex(&regex, regexp_text.c_str());

    // Open directory
    vector<string> files; // list of files in directory
    DIR *dp; // directory pointer
    struct dirent *dirp;
    if((dp  = opendir(blastdir.c_str())) == NULL) {
        cout << "Error: opening cirectory " << blastdir << endl;
        return -1;
    }

    // Add files in directory to files vector
    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);

    // Check to see if files match parsed.txt pattern
    string dbname, tempdbname;
    string dbtag;
    bool usedb ;
    DBTYPE dbtype;
    string (*idextractor)(const char*);

    for (unsigned int i = 0; i < files.size(); i++) {
        dbname = getpattern(&regex, files[i].c_str(), 1) ;
        tempdbname = to_upper(dbname);
        usedb = false;
        if( dbname.size() > 0)  {
            // Pattern found: add database name, filename, and weight to db_info
            // Check for particular databases to add relevant idextractors
            if(tempdbname.find("KEGG") != std::string::npos ) {
               usedb = true;
               dbtype = KEGG;
               idextractor = getKEGGID;
            }

            if(tempdbname.find("COG") != std::string::npos ) {
               usedb = true;
               dbtype = COG;
               idextractor = getCOGID;
            }

            if(tempdbname.find("SEED") != std::string::npos ) {
               usedb = true;
               dbtype = SEED;
               idextractor = getSEEDID;
            }
            
            if(tempdbname.find("CAZY") != std::string::npos ) {
               usedb = true;
               dbtype = CAZY;
               idextractor = getCAZYID;
            }

            if(usedb) {
                db_info.dbtypes.push_back(dbtype); // TODO: May not be unnessisary
                db_info.idextractors[dbname] = idextractor;
            }
            db_info.db_names.push_back(dbname);
            db_info.input_blastouts.push_back(files[i]);
            db_info.weight_dbs.push_back(1.0);
        }
    }

    return 0;
}

/*
 * Extracts ORF_ID from GFF line. Used for doing an on-disk sort.
 */
string orf_extractor_from_gff(const string &line){
   char buf[10000];
   string orfid;
   try{
      orfid  = ShortenORFId(split_n_pick(split_n_pick(split_n_pick(line, buf, '\t', 8), buf, ';',0), buf,'=',1));
    }
    catch(...) {
       return line;
    }
   return orfid;
}

/*
 * Main thread function to create and launch threads for ORF annotation and count reduction
 */
void createThreadsAnnotate(int num_threads, THREAD_DATA_ANNOT *thread_data, WRITER_DATA_ANNOT *writer_data) {
    // Create threads
    pthread_t *threads;
    if( (threads = (pthread_t *) malloc( sizeof(pthread_t) *num_threads) ) == 0 ) {
      std::cout << "Error in allocating threads" << std::endl;
    }
    int rc;
    
    // Launch annotateOrfsForDBs
    for(int i = 0; i < num_threads; i++) {
       if( ( rc = pthread_create(&threads[i], NULL, annotateOrfsForDBs, (void *) (thread_data+i)) ) ) {
         cout << "Error: unable to create thread, " << rc << endl;
         exit(-1);
       }   
    }
    
    cout << "Finished annotateOrfsForDBs()" << endl;
    
    // Join threads
    void *status;
    for(int i=0; i<num_threads; i++) {
      rc = pthread_join(threads[i], &status);
      if (rc) {
         printf("ERROR: return code from pthread_join() is %d\n", rc);
         exit(-1);
      }
    }
    
    cout << "Joined threads from annotateOrfsForDBs()" << endl;
    
    // Launch reduceAnnotations with annotated writer_data
    pthread_t reducer_thread; //(pthread_t *)malloc(sizeof(pthread_t));
    if( ( rc = pthread_create(&reducer_thread, NULL, reduceAnnotations, (void *) (writer_data)) ) ) {
         cout << "Error: unable to create thread," << rc << endl;
         exit(-1);
    }
    
    cout << "Finished reduceAnnotations()" << endl;

    // Join threads
    rc = pthread_join(reducer_thread, &status);
    if (rc) {
         printf("ERROR: return code from pthread_join() is %d\n", rc);
         exit(-1);
    }
    
    cout << "Joined threads from reduceAnnotations()" << endl;
}

/*
 * Given a list of annotations in annotation_list returns the best annotation object with the largest bit score
 */
ANNOTATION* getBestAnnotation(vector<ANNOTATION *> annotation_list, ANNOTATION* annotation) {
    if (annotation_list.size() > 0) {
        annotation = annotation_list[0];
        for(vector<ANNOTATION*>::iterator itr = annotation_list.begin(); itr != annotation_list.end(); ++itr) {
            if ((*itr)->bitscore > annotation->bitscore) {
                annotation = *itr;
            }
        }
    }
    return annotation;
}

/*
 * Main threaded annotation function. Scans ORF annotations and adds their count to for the hierarchical
 * annotation tables
 */
void* annotateOrfsForDBs(void *_data) {

    THREAD_DATA_ANNOT *data = static_cast<THREAD_DATA_ANNOT *> (_data);
    
    // Annotation fields
    ANNOTATION *annotation = new ANNOTATION();
    vector<ANNOTATION *> annotation_list;
    ANNOTATION *metacyc_annotation = new ANNOTATION(); // MetaCyc annotation for matches from MetaCyc Trie
    string (*idextractor) (const char *); // Specific annotation ID extractor
    IDTREE *idtree = new IDTREE(); // Generic DB parser
    unsigned int score; // annotation information score
    string db_id = ""; // database
    string db_name = ""; // helper string for database name
    
    // Ptools Annotation variables
    PTOOLS_NODE *ptools_ptr = NULL;
    bool complete = false;
    vector <string> word_list;
    vector <string> max_word_list;
    string annotation_str = "";
    char buf[10000];
    vector <char *> annotation_words; // Vector of annotation words
    string annotation_product;
    string metacyc_annotation_product;
    bool found_ptools_annotation = false; // Variable to flag of RefSeq annotation (taxa info for later WTD calculation)

    // LCA calculation variables
    NCBITree *ncbi_tree = data->ncbi_tree; // NCBI Taxonomy Tree
    vector<string> ncbi_taxa_ids; // Vector of parsed taxonomy ids
    string db_name_upper = "";
    vector <char *> ncbiID_parsed_fields;
    string lca = "";
    string refseq_id = "";
    char buf2[1000];

    // Main loop
    for( vector<string>::iterator it = data->orfids.begin(); it != data->orfids.end(); it++ )  {
        // For each orf_id

        // Clear previous annotation data
        *annotation = ANNOTATION(); // clear
        db_name = "";

        for( unsigned int j = 0; j < data->db_info.db_names.size(); j++ ) { 
            // For each database j
            db_name = data->db_info.db_names[j];
            db_name_upper = to_upper(db_name);
            found_ptools_annotation = false;
            *metacyc_annotation = ANNOTATION();

            // Add database to map if new
            if (data->dbNamesToHierachyIdentifierCounts.find(db_name) == data->dbNamesToHierachyIdentifierCounts.end()) {
                map<string, int> newHierarchyMap;
                data->dbNamesToHierachyIdentifierCounts[db_name] = newHierarchyMap;
            }

            // Get annotation if orf_id found
            if( data->annot_objects[db_name].find(*it) != data->annot_objects[db_name].end()) {
                // Annotation orf_id found in database j
                annotation = new ANNOTATION();

                // Get the hit with the highest score
                annotation_list = data->annot_objects[db_name][*it];
                annotation = getBestAnnotation(annotation_list, annotation);
                annotation->dbname = db_name;

                // Pass annotation through ptools tree and add if hit
                annotation_product = annotation->product;

                // Split annotation into separate words
                split(annotation_product, annotation_words, buf, ' ');

                // Process annotation through MetaCyc trie
                metacyc_annotation_product = processAnnotationForPtools(annotation_words,
                                                                        data->root,
                                                                        ptools_ptr,
                                                                        complete,
                                                                        word_list,
                                                                        max_word_list,
                                                                        annotation_str);

                // Calculate information annotation score
                score = computeAnnotationValue(annotation) * data->db_info.weight_dbs[j];
                annotation->annotation_score = score;

                // If MetaCyc or EC number found add annotation to metaCycHits for ptools
                if (metacyc_annotation_product != "") {
                    // Matched some MetaCyc annotation
                    *metacyc_annotation = *annotation;
                    metacyc_annotation->product = metacyc_annotation_product;
                    metacyc_annotation->ptools_match = true;
                    data->metaCycHits.push_back(*metacyc_annotation);
                    found_ptools_annotation = true;
                } else if ( annotation->ec != "" ) {
                    // Valid EC number
                    *metacyc_annotation = *annotation;
                    annotation->ptools_match = false;
                    data->metaCycHits.push_back(*annotation);
                    found_ptools_annotation = true;
                }
                
                // Extract db_id from annotation text
                if (data->db_info.idextractors.find(db_name) != data->db_info.idextractors.end()) {
                    // Use specific idextractor function if found
                    idextractor = data->db_info.idextractors[db_name];
                    db_id = idextractor(annotation->product.c_str());
                    if (db_id != "") {
                        if (data->dbNamesToHierachyIdentifierCounts[db_name].find(db_id) == data->dbNamesToHierachyIdentifierCounts[db_name].end()) {
                            // First time seeing db
                            data->dbNamesToHierachyIdentifierCounts[db_name][db_id] = 0;
                        }
                        data->dbNamesToHierachyIdentifierCounts[db_name][db_id]++; // add count to identifier
                    }
                } else if ( (data->dbNamesToHierarchyIdTree.find(db_name) != data->dbNamesToHierarchyIdTree.end())) {
                    // Generic id extractor case
                    idtree = data->dbNamesToHierarchyIdTree[db_name];
                    if (idtree->find(annotation->product).size() > 0) {
                        db_id = idtree->find(annotation->product);
                        if (db_id != "") {
                            if (data->dbNamesToHierachyIdentifierCounts[db_name].find(db_id) == data->dbNamesToHierachyIdentifierCounts[db_name].end()) {
                                // First time seeing db
                                data->dbNamesToHierachyIdentifierCounts[db_name][db_id] = 0;
                            }
                            
                            data->dbNamesToHierachyIdentifierCounts[db_name][db_id]++;
                        }
                    }
                } else if (db_name_upper.find("REFSEQ") != std::string::npos) {
                    // RefSeq case

                    // Calculate LCA
                    ncbi_taxa_ids.clear(); // Collection of NCBI IDs
                    for (vector<ANNOTATION*>::iterator itr = annotation_list.begin(); itr != annotation_list.end(); ++itr) {
                        ncbiID_parsed_fields.clear();
                        refseq_id = (*itr)->target;
                        split(refseq_id, ncbiID_parsed_fields, buf2, '|');
                        refseq_id = ncbiID_parsed_fields[3];
                        if (ncbi_tree->RefSeqID_to_NCBI_ID.find(refseq_id) != ncbi_tree->RefSeqID_to_NCBI_ID.end()) {
                            ncbi_taxa_ids.push_back(ncbi_tree->RefSeqID_to_NCBI_ID[refseq_id]);
                        }
                    }

                    lca = ncbi_tree->getLCA(ncbi_taxa_ids); // search tree

                    // Associate LCA to ptools annotation (refseq only)
                    if (found_ptools_annotation) {
                        // Get ptools annotaiton string
                        string ptools_key = "";
                        if(metacyc_annotation->ptools_match) {
                            ptools_key = metacyc_annotation->product;
                        } else {
                            ptools_key = metacyc_annotation->ec;
                        }

                        if (ptools_key != "") {
                            // Add string to map if needed
                            if (data->metaCycHitToNCBITaxonomy.find(ptools_key) == data->metaCycHitToNCBITaxonomy.end()) {
                                data->metaCycHitToNCBITaxonomy[ptools_key] = map<string, int>();
                            }

                            if (data->metaCycHitToNCBITaxonomy[ptools_key].find(lca) == data->metaCycHitToNCBITaxonomy[ptools_key].end()) {
                                data->metaCycHitToNCBITaxonomy[ptools_key][lca] = 0;
                            }
                            data->metaCycHitToNCBITaxonomy[ptools_key][lca]++;
                        }
                    }

                    // Add LCA to the Hierarchy count
                    if (data->dbNamesToHierachyIdentifierCounts[db_name].find(lca) == data->dbNamesToHierachyIdentifierCounts[db_name].end()) {
                        // First time seeing db_id
                        data->dbNamesToHierachyIdentifierCounts[db_name][lca] = 0;
                    }
                    data->dbNamesToHierachyIdentifierCounts[db_name][lca]++; // Add count to identifier
                }
            }
        }
    }

    return (void *) NULL;
}


/*
 * Collects annotation counts from each thread and combines them into their respective global data structures in writer_data
 * to be written out to disk in the main function. Runs after each batch of annotations summarized.
 * dbNamesToHierachyIdentifierCounts: database_name -> db_id -> count
 * globalMetaCycNamesToAnnotations: <MetaCycProductName/ECNumber> -> ANNOTATION object
 * globalMetaCycNamesToDbCounts: <MetaCycProductName/ECNumber> -> database -> count
 */
void *reduceAnnotations( void *_writer_data) {

    WRITER_DATA_ANNOT *writer_data = (WRITER_DATA_ANNOT *) _writer_data;
    THREAD_DATA_ANNOT *thread_data = writer_data->thread_data; // get annotation data
    unsigned int num_threads = writer_data->num_threads;
    
    for(unsigned int i = 0; i < num_threads; i++) {
        // for each thread
        
        // reduce DB hierarchy counts in dbNamesToHierachyIdentifierCounts
        std::cout << "Reducing DbNamesToHierachyIdentifierCounts results from thread " << i << std::endl;
        for ( map<string, map<string, int> >::iterator db_itr = thread_data[i].dbNamesToHierachyIdentifierCounts.begin();
              db_itr != thread_data[i].dbNamesToHierachyIdentifierCounts.end();
              db_itr++
              ) {

            // For each db_name
            if ( writer_data->globalDbNamesToHierachyIdentifierCounts.find(db_itr->first) == writer_data->globalDbNamesToHierachyIdentifierCounts.end()) {
                // create database in map if not present
                map<string, int> newHierarchyMap;
                writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first] = newHierarchyMap;
            }

            // Iterate through all ids
            for (map<string, int>::iterator id_itr = thread_data[i].dbNamesToHierachyIdentifierCounts[db_itr->first].begin();
                 id_itr != thread_data[i].dbNamesToHierachyIdentifierCounts[db_itr->first].end();
                 id_itr++) {
                 if ( writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first].find(id_itr->first) == writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first].end()) {
                     // create id slot if not present in database
                     writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first][id_itr->first] = 0;
                 }
                 // add to global count
                 writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first][id_itr->first] += id_itr->second;
            }
        }

        // Calculate the global MetaCyc annotations
        string metacyc_label = ""; // EC number match or MetaCyc annotation used as key

        // For each annotation
        for (vector<ANNOTATION>::iterator mc_itr = thread_data[i].metaCycHits.begin();
             mc_itr != thread_data[i].metaCycHits.end();
             mc_itr++) {
             if (mc_itr->ptools_match) {
                 metacyc_label = mc_itr->product;

             } else {
                 metacyc_label = mc_itr->ec;
             }
             // Create new entrys if needed
             // Set annotation (default is the first)
             if (writer_data->globalMetaCycNamesToAnnotations.find(metacyc_label) == writer_data->globalMetaCycNamesToAnnotations.end()) {
                 writer_data->globalMetaCycNamesToAnnotations[metacyc_label] = *mc_itr;
             }
             if (writer_data->globalMetaCycNamesToDbCounts.find(metacyc_label) == writer_data->globalMetaCycNamesToDbCounts.end()) {
                 map<string, int> newMetaCycDbNameCountMap;
                 writer_data->globalMetaCycNamesToDbCounts[metacyc_label] = newMetaCycDbNameCountMap;
             }

             // Add database count
             if (writer_data->globalMetaCycNamesToDbCounts[metacyc_label].find(mc_itr->dbname) == writer_data->globalMetaCycNamesToDbCounts[metacyc_label].end()) {
                 writer_data->globalMetaCycNamesToDbCounts[metacyc_label][mc_itr->dbname] = 0;
             }
             writer_data->globalMetaCycNamesToDbCounts[metacyc_label][mc_itr->dbname]++;
        }

        // Calculate global NCBI Taxonomy counts
        metacyc_label = "";
        string taxa_string = "";
        for (map< string, map< string, int > >::iterator itr = thread_data[i].metaCycHitToNCBITaxonomy.begin();
                itr != thread_data[i].metaCycHitToNCBITaxonomy.end();
                itr++) {
            metacyc_label = itr->first;
            if (writer_data->globalMetaCycHitToNCBITaxonomy.find(metacyc_label) ==
                writer_data->globalMetaCycHitToNCBITaxonomy.end()) {
                writer_data->globalMetaCycHitToNCBITaxonomy[metacyc_label] = map<string, int>();
            }

            for (map<string, int>::iterator itr2 = itr->second.begin();
                 itr2 != itr->second.end();
                 itr2++) {
                taxa_string = itr2->first;
                if (writer_data->globalMetaCycHitToNCBITaxonomy[metacyc_label].find(taxa_string) ==
                    writer_data->globalMetaCycHitToNCBITaxonomy[metacyc_label].end()) {
                    writer_data->globalMetaCycHitToNCBITaxonomy[metacyc_label][taxa_string] = 0;
                }
                writer_data->globalMetaCycHitToNCBITaxonomy[metacyc_label][taxa_string] += itr2->second;

            }
        }
        thread_data[i].clear();
    }
    
    return (void *)NULL;

}

/*
 * Print out MetaCyc word trie using a stack-based recursive pre-order traversal
 */
void printMetaCycTree(PTOOLS_NODE *root) {
    vector<PTOOLS_NODE *> node_stack; // Stack of nodes
    vector<string> path; // Stack for path tracing
    string headpath = ""; // Path string

    // Push first nodes
    node_stack.push_back(root);
    path.push_back(headpath);

    while(!node_stack.empty()) {
        // Pop stacks
        PTOOLS_NODE* top = node_stack.back();
        node_stack.pop_back();
        headpath = path.back(); // Get current path
        path.pop_back();

        if (top->complete) {
            // Complete MetaCyc annotation
            cout << headpath << endl;
        }
        // Push children to stacks
        map<string, PTOOLS_NODE*> children = top->children;
        for (std::map<string, PTOOLS_NODE*>::iterator itr=children.begin(); itr != children.end(); ++itr) {
            node_stack.push_back(itr->second);
            path.push_back(headpath + itr->first + " "); // Keep track of path
        }
    }
}

/*
 * Scans Pathway Tools reaction file and puts reaction strings into ptools_list
 */
void processPToolsRxnsFile( string ptools_rxn_file, vector<string> &ptools_list ) {
    string filename  = ptools_rxn_file;
    ifstream input;
    std::cout << "Reading ptools_rxn_file " << filename <<  std::endl;
    input.open(filename.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
        return ;
    }

    // Scanning variables
    int count =0;
    string line;
    char buf[10000]; // Temp buffer
    vector <char *> fields; // Vector for parsed fields

    while( std::getline(input, line ).good()) {
        // Skip header
        if (count == 0) {
            count++;
            continue;
        }

        split(line, fields, buf,'\t');
        ptools_list.push_back(fields[1]); // Ptools enzyme description field.
        count++;
    }

    input.close();
}

/*
 * Checks an individual annotation (list of words) and its word-strings for 'complete' status in the MetaCyc trie. If a
 * 'complete' annotation node reached.
 */
string processAnnotationForPtools(vector <char *> annotation_words, PTOOLS_NODE *root, PTOOLS_NODE *ptools_ptr,
                                  bool complete, vector <string> word_list, vector <string> max_word_list,
                                  string annotation) {
    complete = false;
    word_list.clear();
    max_word_list.clear();
    annotation.clear();

    // Try to push current word
    for (unsigned int i = 0; i < annotation_words.size(); i++) {
        word_list.clear();
        ptools_ptr = root; // reset to root
        for (unsigned int j = i; j < annotation_words.size(); j++) {
            string word = string(annotation_words[j]);
            if (pushWordForward(word, ptools_ptr)) {
                // Check to see if pointer now at word that completes an annotation
                word_list.push_back(ptools_ptr->id2);
                if(ptools_ptr->complete) {
                    complete = true;
                    if (word_list.size() > max_word_list.size()) {
                        max_word_list = word_list;
                    }
                }
            } else {
                //cout << "Not Found " << annotation_words[j] << endl;
                break;
            }
        }
    }
    if (complete) {
        // Found complete ptools annotation in annotation_words
        vector <string>::iterator itr;
        for(itr = max_word_list.begin(); itr != max_word_list.end(); itr++) {
            annotation = annotation + *itr + " ";
        }

        return annotation.substr(0, annotation.size()-1);
        // cout << "Found complete" << endl;
        // cout << "Annotation: " << annotation << endl;
    }
    return "";

}

/*
 * Checks to see if current word is a child of the current node pointed at by ptools_ptr. Moves pointer to child and
 * returns true if so, and returns false otherwise.
 */
bool pushWordForward(string word, PTOOLS_NODE *& ptools_ptr) {
    if( ptools_ptr->hasChild(word) ) {
        // Word is a child
        ptools_ptr = ptools_ptr->getChildNode(word);
        return true;
    } else {
        return false;
    }
}

/*
 * Function takes all ingredients for a Pathway Tools pf file entry, creates the string, and writes the string to the
 * .pf output file.
 */
void writePfEntry(string orf_id, string annotation_product, string ec_number, string start_base, string end_base, std::ofstream &output) {

    // Example .pf format
    // ID <orf_id>
    // NAME <orf_id>
    // STARTBASE <start_base>
    // ENDBASE <start_base>+<length>
    // PRODUCT-TYPE P
    // FUNCTION <annotation_product>
    // EC <ec_number>
    // METACYC <rxn_id> (optional)
    // PRODUCT-ID <product_id> (optional)
    // //

    // remove last string
    // cout << "start_base: " << start_base << " end_base: " << end_base << endl;

    string pf_entry = "";
    pf_entry = pf_entry + "ID\t" + orf_id + "\n";
    pf_entry = pf_entry + "NAME\t" + orf_id + "\n";
    pf_entry = pf_entry + "STARTBASE\t" + start_base + "\n";
    pf_entry = pf_entry + "ENDBASE\t" + end_base + "\n";
    pf_entry = pf_entry + "PRODUCT-TYPE\tP" + "\n";
    pf_entry = pf_entry + "FUNCTION\t" + annotation_product + "\n"; // remove extract character
    if (ec_number != "")
        pf_entry = pf_entry + "EC\t" + ec_number + "\n";
    pf_entry = pf_entry + "//" + "\n";
    output << pf_entry;
}


/*
 * Creates string of taxonomies and counts for output
 */
string getRefSeqTaxonomiesForPtools(string &metacyc_anno,
                                    WRITER_DATA_ANNOT* writer_data,
                                    string &taxaline,
                                    string &full_name,
                                    string &taxa) {
    taxaline = "";
    full_name = "";
    taxa = "";
    if (writer_data->globalMetaCycHitToNCBITaxonomy.find(metacyc_anno) != writer_data->globalMetaCycHitToNCBITaxonomy.end()) {

        map<string, int> taxacounts = writer_data->globalMetaCycHitToNCBITaxonomy[metacyc_anno];
        taxaline = "{ ";
        int lin_count = 0;
        for (map<string, int>::iterator itr = taxacounts.begin(); itr != taxacounts.end(); itr++) {
            taxa = itr->first;
            ostringstream ss3;
            ss3 << itr->second;
            full_name = "";
            if (writer_data->ncbi_tree->NCBI_ID_to_Common.find(taxa) != writer_data->ncbi_tree->NCBI_ID_to_Common.end()) {
                full_name = writer_data->ncbi_tree->NCBI_ID_to_Common[taxa];
            }
            taxa = full_name + " (" + taxa + ")";
            if (lin_count == 0) {
                taxaline = taxaline + taxa + ":" + ss3.str();
            } else {
                taxaline = taxaline + ", " + taxa + ":" + ss3.str();
            }
            lin_count++;
            ss3.clear();
        }

        taxaline = taxaline + " }";
    }

    return taxaline;
}


/*
 * Creates the inputs for the ptools folder from the writer_data:
 * 0.pf file, ptools_counts.txt, genetic-elements.dat, organism-params.dat
 */
void writePToolsResults(WRITER_DATA_ANNOT* writer_data, string ptools_dir, string sample_name) {

    // Create 0.pf file
    string pf_file = writer_data->options.ptools_dir + "/" + "0.pf";
    ofstream pf_output;
    pf_output.open(pf_file.c_str(), std::ofstream::out);
    if(!pf_output.good()){
        cerr << "Error opening '" << pf_file << "'. " << endl;
        exit(-1);
    }
    int start_base = 0;
    int end_base = 0;
    char start_base_str[30];
    char end_base_str[30];
    string derived_orf_id;

//    map<string, string> anno_to_frameid = map<string, string>();

    // Writeout metaCycHits to 0.pf file
    int i = 0;
    for (map<string, ANNOTATION>::iterator mc_itr =  writer_data->globalMetaCycNamesToAnnotations.begin();
         mc_itr != writer_data->globalMetaCycNamesToAnnotations.end();
         mc_itr++) {
         ANNOTATION anno = mc_itr->second;
         end_base = start_base + anno.length;

         sprintf(start_base_str, "%d", start_base);
         sprintf(end_base_str, "%d", end_base);
         ostringstream ss;
         ss << i;
         derived_orf_id = "DIR_" + sample_name + "_" + ss.str(); // modified orf_id

//         if (anno_to_frameid.find(derived_orf_id) == anno_to_frameid.end()) {
//             anno_to_frameid[derived_orf_id] = anno.product;
//         }

         writePfEntry(derived_orf_id,
                      anno.product,
                      anno.ec,
                     string(start_base_str),
                     string(end_base_str),
                     pf_output);
         // update startbase
         start_base = start_base + anno.length + 10;
         i++;
         ss.clear();
    }

    pf_output.close();

    // Create counts file
    string counts_file = writer_data->options.ptools_dir + "/" + "ptools_counts.txt";

    pf_output.open(counts_file.c_str(), std::ofstream::out);
    if(!pf_output.good()){
        cerr << "Error opening '" << pf_file << "'. " << endl;
        exit(-1);
    }

    // Header line with database names
    string header_line = "Annotation";
    for (unsigned int i = 0; i < writer_data->db_info.db_names.size(); ++i) {
        header_line = header_line + "\t" + writer_data->db_info.db_names[i];
    }
    header_line = header_line + "\t" + "Total" + "\t" + "Taxa";
    header_line = header_line + "\n";
    pf_output << header_line;

    string taxaline = "";
    string full_name = "";
    string taxa = "";
    string metacyc_anno = "";

    // Create count line
    for (map<string, map<string, int> >::iterator mc_itr =  writer_data->globalMetaCycNamesToDbCounts.begin();
         mc_itr != writer_data->globalMetaCycNamesToDbCounts.end();
         mc_itr++) {
         string line = "";
         metacyc_anno = "";
         line = mc_itr->first;
         metacyc_anno = mc_itr->first;

         int total = 0;
         map<string, int> db_counts = mc_itr->second;
         for (unsigned int i = 0; i < writer_data->db_info.db_names.size(); ++i) {
             if (db_counts.find(writer_data->db_info.db_names[i]) == db_counts.end() ) {
                 line = line + "\t" + "0";
             } else {
                 ostringstream ss1;
                 ss1 << db_counts[writer_data->db_info.db_names[i]];
                 line = line + "\t" + ss1.str();
                 total += db_counts[writer_data->db_info.db_names[i]];
             }
         }

         string taxonomies = getRefSeqTaxonomiesForPtools(metacyc_anno, writer_data, taxaline, full_name, taxa);
         ostringstream ss2;
         ss2 << total;
         line = line + "\t" + ss2.str();
         line = line + "\t" + taxonomies;
         line = line + "\n";
         pf_output << line;
    }
    pf_output.close();

    // Create genetic-elements.dat
    // ID      0
    // NAME    0
    // TYPE    :READ/CONTIG
    // ANNOT-FILE      0.pf

    string genetic_elements_file = ptools_dir + "/" + "genetic-elements.dat";
    ofstream output;
    output.open(genetic_elements_file.c_str(), std::ofstream::out);
    if(!output.good()){
        cerr << "Error opening '" << genetic_elements_file << "'. " << endl;
        return ;
    }

    string genetic_elements_text = "";
    genetic_elements_text = genetic_elements_text + "ID\t0" + "\n";
    genetic_elements_text = genetic_elements_text + "NAME\t0" + "\n";
    genetic_elements_text = genetic_elements_text + "TYPE\t:READ/CONTIG" + "\n";
    genetic_elements_text = genetic_elements_text + "ANNOT-FILE\t0.pf" + "\n";

    output << genetic_elements_text;

    output.close();

    // Create organism-params.dat
    // ID      hmp_airways_SRS014682
    // STORAGE FILE
    // NAME    hmp_airways_SRS014682
    // ABBREV-NAME     hmp_airways_SRS014682
    // STRAIN  1
    // RANK    |species|
    // NCBI-TAXON-ID   12908

    string orgamism_params_file = ptools_dir + "/" + "organism-params.dat";
    output.open(orgamism_params_file.c_str(), std::ofstream::out);
    if(!output.good()){
        cerr << "Error opening '" << orgamism_params_file << "'. " << endl;
        return ;
    }

    string organism_params_text = "";
    organism_params_text = organism_params_text + "ID\t" + sample_name + "\n";
    organism_params_text = organism_params_text + "STORAGE FILE" + "\n";
    organism_params_text = organism_params_text + "NAME\t" + sample_name + "\n";
    organism_params_text = organism_params_text + "ABBREV-NAME\t" + sample_name + "\n";
    organism_params_text = organism_params_text + "STRAIN\t1" + "\n";
    organism_params_text = organism_params_text + "RANK\t|species|" + "\n";
    organism_params_text = organism_params_text + "NCBI-TAXON-ID\t12908" + "\n";

    output << organism_params_text;

    output.close();
}
 
bool createFunctionWeights(const string &inputfile, const string &outputfile) {

    std::cout << "file to use " <<  inputfile << std::endl;
    std::ifstream input;
    // char buf[1000];
 
    input.open(inputfile.c_str(), std::ifstream::in);
    if(!input.good()){
       std::cerr << "Error opening '"<< inputfile <<"'. Bailing out." << std::endl;
       return false;
    }   
      // int *counts  = (int *)calloc(options.num_threads, sizeof(int));
    string function, line;
    string prevfunc="";
    bool isfirst = true;
    unsigned int count = 1;


    std::ofstream output;
    output.open(outputfile.c_str(), std::ifstream::out);

    while( std::getline(input, line ).good()) {
        function = line;
         
        if( prevfunc == function)    
             count++;
        else
           if(isfirst ) {
              isfirst = false;
           }
           else {
             output  << prevfunc << "\t" << count << std::endl;
             count =1;
           }
        
        prevfunc = function;
     }   

    input.close();
    output.close();

    return true;

}

