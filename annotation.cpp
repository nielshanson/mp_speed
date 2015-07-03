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

int processParsedBlastout(string db_name, float weight, string blastoutput, MPAnnotateOptions options, map<string, ANNOTATION> &annotation_results) {

    if (options.debug) cout << "In processParsedBlastout()" << endl;

    // Prepare inputs and buffers
    string filename = options.blast_dir + "/" + blastoutput; // BLAST/LASTout.parsed.txt
    std::ifstream input; // Input filestream
    char buf[10000]; // Temp buffer TODO check to see if multiple buffers nessisary
    vector <char *> fields; // Vector for parsed fields

    int count = 0; // Line count

    // char tempbuf[1000];

    if (options.debug) {
        cout << "Reading ParsedBlastout: " << filename << "\n";
    }

    // Open the file.
    input.open(filename.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '" << filename << "'. Bailing out." << std::endl;
        exit(1);
    }

    // ANNOTATION fields for parsing BLAST/LASTout.annotated.txt files
    ANNOTATION annotation;
    string line;
    string query_id;
    string bsr;
    string ec;
    string product;
    string taxonomy = "";

    // Parse annotation hit line, create ANNOTATION objects
    while( std::getline(input, line ).good()) {
        split(line, fields, buf,'\t');
        if (count == 0) { count++; continue; };

        if( fields.size()!=10) {
            // Not a parsed blast file
            cerr << "Parsed BLAST/LASTout file " <<  filename << " did not have the 10 columns." << endl;
            query_id.clear();
            input.close();
            return 1;
        }

        query_id = fields[0];
        bsr = fields[4];
        ec = fields[8];
        product = fields[9];

        db_name = to_upper(db_name);

        // Construct annotation
        ANNOTATION annotation;
        annotation.bsr = atof(bsr.c_str());
        annotation.ec = ec;
        annotation.product = product;

        // Extract RefSeq taxonomy from product field
        if (options.taxonomy && (db_name.find("REFSEQ") != std::string::npos)) {
            // TODO: Could be more reliable to get taxonomy from GI number.
            taxonomy = getTaxonomyFromProduct(product.c_str());
        }
        annotation.taxonomy = taxonomy;

        // Compute information score of current annotation.
        annotation.value = computeAnnotationValue(&annotation) * weight;

        // Check to see if current query_id is already in this database's annotation_results
        if (annotation_results.count(query_id) <= 0) {
            annotation_results[query_id] = annotation;
        }
        else {
            // Replace existing ANNOTATION with the current annotation if strong annotation score
            if (annotation_results[query_id].value < annotation.value) {
                annotation_results[query_id] = annotation;
            }
        }

        if (options.debug) {
            if (count%PRINT_INTERVAL==0)
                std::cout << "x " << count << std::endl;
        }

        count++;
    }
    input.close();

    if (options.debug) {
        std::cout << "Number of parsed BLAST/LAST results loaded " <<  count << std::endl;
    }
    return count;
}


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
                db_info.db_names.push_back(dbname);
                db_info.input_blastouts.push_back(files[i]);
                db_info.weight_dbs.push_back(1.0);
                db_info.dbtypes.push_back(dbtype); // TODO: May not be unnessisary
                db_info.idextractors[dbname] = idextractor;
            } else {
                db_info.db_names.push_back(dbname);
                db_info.input_blastouts.push_back(files[i]);
                db_info.weight_dbs.push_back(1.0);
            }
        }
    }

    return 0;
}

void createAnnotation(map<string, float> dbname_weight, ANNOTATION_RESULTS results_dictionary, MPAnnotateOptions options, map<string, unsigned int> contig_lengths) {
    // create_annotation(dbname_weight, results_dictionary, opts.input_gff, opts.rRNA_16S, opts.tRNA, opts.output_gff, opts.output_comparative_annotation, contig_lengths)
    // orf_dictionary={};

    if (options.debug) {
        cout << "In createAnnotation()" << endl;
    }
    // string input_gff = options.input_gff;
    // string rRNA_16S = options.rRNA_16S;
    // string tRNA = options.tRNA;
    // string output_gff = options.output_gff;
    // string output_comp_annot = options.output_comp_annot;

    // // read input ORFs the input GFF file
    // cout << options.input_gff << endl;

    // // Prepare inputs and buffers
    // string filename = options.blast_dir + "/" + input_gff; // BLAST/LASTout.parsed.txt

    // std::ifstream input; // Input filestream
    // vector <char *> fields; // Vector for parsed fields





//    output_gff_tmp = output_gff + ".tmp";
//    outputgff_file = open( output_gff_tmp, 'w');
//    output_comp_annot_file1 = open( output_comparative_annotation + '.1.txt', 'w');
//    output_comp_annot_file2 = open( output_comparative_annotation + '.2.txt', 'w');
//
//    output_comp_annot_file1_Str = 'orf_id\tref dbname\tEC\tproduct\tvalue';
//    fprintf(output_comp_annot_file1,'%s\n', output_comp_annot_file1_Str);
//
//    output_comp_annot_file2_Str = 'orf_id'
//    dbnames = dbname_weight.keys()
//    for dbname in dbnames:
//    weight = dbname_weight[dbname]
//    output_comp_annot_file2_Str += '\t{0}(EC) \t{0}(product)\t{0}(value)'.format(dbname)
//    fprintf(output_comp_annot_file2,'%s\n', output_comp_annot_file2_Str)
}


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
 * Create threads and run annotateOrfsForDBs and writeAnnotatedGFFs 
 */
void createThreadsAnnotate(int num_threads, THREAD_DATA_ANNOT *thread_data, WRITER_DATA_ANNOT *writer_data) {
    // Create threads
    pthread_t *threads;
    if( (threads = (pthread_t *) malloc( sizeof(pthread_t) *num_threads) ) == 0 ) {
      std::cout << "Error in allocating threads" << std::endl;
    }
    int rc;
    
    // Lanuch annotateOrfsForDBs
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
    
    cout << "Joined threads..." << endl;
    
    // Create new writer_tread for writeannotate GFFs
    pthread_t writer_thread; //(pthread_t *)malloc(sizeof(pthread_t)); 
    if( ( rc = pthread_create(&writer_thread, NULL, writeAnnotatedGFFs, (void *) (writer_data)) ) ) {
         cout << "Error: unable to create thread," << rc << endl;
         exit(-1);
    }
    
    cout << "Finished annotateOrfsForDBs()" << endl;

    rc = pthread_join(writer_thread, &status);
    if (rc) {
         printf("ERROR: return code from pthread_join() is %d\n", rc);
         exit(-1);
    }
    
    cout << "Joined threads..." << endl;
    
}

/*
 * Creates annotation for functional_and_taxonomic table using hits from available
 * databases.
 */
void *annotateOrfsForDBs( void *_data) {
    // cast _data to THREAD_DATA_ANNOT
    THREAD_DATA_ANNOT *data = static_cast<THREAD_DATA_ANNOT *> (_data);
    
    // annotation fields
    unsigned int max_score = 0, score;
    string db_id = "";
    
    // process fields
    bool success;
    ANNOTATION *annotation = new ANNOTATION();
    ANNOTATION *final_annotation = new ANNOTATION();
    string (*idextractor) (const char *);
    vector<string>::iterator it;
    unsigned int total_count = 0;
    unsigned int annotated_count = 0;
    IDTREE *idtree = new IDTREE();
    
    // Ptools Annotation variables
    PTOOLS_NODE *ptools_ptr; 
    bool complete = false;
    vector <string> word_list;
    vector <string> max_word_list;
    string annotation_str = "";
    char buf[10000]; // Temp buffer
    vector <char *> annotation_words; // Vector of annotation words
    string annotation_product;
    string metacyc_annotation_product;
    
    for( it = data->orfids.begin(); it != data->orfids.end(); it++ )  {
        // for each ORF 
        total_count++;
        max_score = 0; // reset score
        success = false;
        *annotation = ANNOTATION(); // clear
        *final_annotation = ANNOTATION();
        
        string db_name = ""; // helper string
        
        for( unsigned int j = 0; j < data->db_info.db_names.size(); j++ ) { 
            // For each database j
            db_name = data->db_info.db_names[j];
            if (data->dbNamesToHierachyIdentifierCounts.find(db_name) == data->dbNamesToHierachyIdentifierCounts.end()) {
                // add database to map
                map<string, int> newHierarchyMap;
                data->dbNamesToHierachyIdentifierCounts[db_name] = newHierarchyMap;
            }
            
            if( data->annot_objects[db_name].find(*it) != data->annot_objects[db_name].end()) {
                // Annotation orf_id found in database j
                success = true;
                // Get the annotation
                annotation = data->annot_objects[db_name][*it];
                score = computeAnnotationValue(annotation) * data->db_info.weight_dbs[j]; // calculate information annotation score
                
                // Set annotation to best hit
                if (score > max_score) {
                    // set final annotation to main fields
                    final_annotation->orf_id = *it;
                    final_annotation->product = annotation->product;
                    final_annotation->bsr = annotation->bsr;
                    final_annotation->value = annotation->value;
                    final_annotation->ec = annotation->ec;
                    final_annotation->taxonomy = annotation->taxonomy;
                    final_annotation->length = annotation->length;
                    
                    max_score = score; // update score
                }
                
                // Get db_id from annotation if need
                if (data->db_info.idextractors.find(db_name) != data->db_info.idextractors.end()) {
                    // Use idextractor function if found
                    idextractor = data->db_info.idextractors[db_name];
                    db_id = idextractor(annotation->product.c_str());
                    if (data->dbNamesToHierachyIdentifierCounts[db_name].find(db_id) == data->dbNamesToHierachyIdentifierCounts[db_name].end()) {
                        // First time seeing db
                        data->dbNamesToHierachyIdentifierCounts[db_name][db_id] = 0;
                    }
                    data->dbNamesToHierachyIdentifierCounts[db_name][db_id] ++; // add count to identifier
                    
                    final_annotation->db_ids[db_name] = db_id; // TODO: probably not needed anymroe
                } else if (data->dbNamesToHierarchyIdTree.find(db_name) != data->dbNamesToHierarchyIdTree.end()) {
                    // Use generic idtree
                    idtree = data->dbNamesToHierarchyIdTree[db_name];
                    if (idtree->find(annotation->product).size() > 0) {
                        db_id = idtree->find(annotation->product);
                        if (data->dbNamesToHierachyIdentifierCounts[db_name].find(db_id) == data->dbNamesToHierachyIdentifierCounts[db_name].end()) {
                            // First time seeing db
                            data->dbNamesToHierachyIdentifierCounts[db_name][db_id] = 0;
                        }
                        
                        data->dbNamesToHierachyIdentifierCounts[db_name][db_id] ++;
                    }
                }
            }
        }
        // TODO: place to check for pathway tools annotation
        if (success) {
            // data->db_hits[*it] = *final_annotation;
            
            // Clear annotation fields
            annotation_product = final_annotation->product;
    
            // Split annotation into separate words
            split(annotation_product, annotation_words, buf, ' ');
    
            // Process annotation through MetaCyc trie
            // annotation_product = processAnnotationForPtools(annotation_words, root);
            metacyc_annotation_product = processAnnotationForPtools(annotation_words, data->root, ptools_ptr, complete, word_list, max_word_list, annotation_str);

            if (metacyc_annotation_product != "") {
                final_annotation->product = metacyc_annotation_product;
                data->metaCycHits.push_back(*final_annotation);
            } else if ( final_annotation->ec != "" ) {
                data->metaCycHits.push_back(*final_annotation);
            }
            // If valid annotation product or EC number
            // if (ec_number != "") {
            //     // If ec_number valid use original annotation
            //     writePfEntry(orf_id, annotation_product, ec_number, start_base, length, output);
            // } else if (metacyc_annotation_product != "") {
            //     writePfEntry(orf_id, metacyc_annotation_product, ec_number, start_base, length, output);
            // }
            
            // metacyc_annotation_product = processAnnotationForPtools(annotation_words, root, ptools_ptr, complete, word_list, max_word_list, annotation);
            // annotated_count++;
        }
    }

    // cout << "annotateOrfsForDBs(): " << annotated_count << " of " << total_count  << " ORFs annotated" << endl;
    return (void *) NULL;

}

//void *writeAnnotatedPreamble(void *_writer_data) {
//
//    char buf[10000];
//    string str;
//    WRITER_DATA_ANNOT *writer_data = (WRITER_DATA_ANNOT *)_writer_data;
//    std::cout << "preamble \n";
//    for(unsigned int i = 0; i < writer_data->db_info.db_names.size(); i++) {
//       writer_data->output[i] <<  "#"<< writer_data->db_info.db_names[i] << std::endl;
//    }
//}

void *writeAnnotatedGFFs( void *_writer_data) {
    
    unsigned int b;  
    // char buf[10000];
    WRITER_DATA_ANNOT *writer_data = (WRITER_DATA_ANNOT *) _writer_data;
    unsigned int num_threads = writer_data->num_threads;
    THREAD_DATA_ANNOT *thread_data = writer_data->thread_data; // get annotation data
    ANNOTATION result_annotation;
    string print_line; // line to print
    
    
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
    
    for(unsigned int i = 0; i < num_threads; i++) {
        // for each thread
        
        b =  (thread_data[i].b+1)%2;
        
        // reduce dbNamesToHierachyIdentifierCounts
        std::cout << "Reducing DbNamesToHierachyIdentifierCounts results from thread " << i << " buffer " << b << std::endl;
        for ( map<string, map<string, int> >::iterator db_itr = thread_data[i].dbNamesToHierachyIdentifierCounts.begin();
              db_itr != thread_data[i].dbNamesToHierachyIdentifierCounts.end();
              db_itr ++
              ) {
            // For each db_name
            if ( writer_data->globalDbNamesToHierachyIdentifierCounts.find(db_itr->first) == writer_data->globalDbNamesToHierachyIdentifierCounts.end()) {
                // create database map if not present in globalDbNamesToHierachyIdentifierCounts
                map<string, int> newHierarchyMap;
                writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first] = newHierarchyMap;
            }
            // Iterate through all ids
            for (map<string, int>::iterator id_itr = thread_data[i].dbNamesToHierachyIdentifierCounts[db_itr->first].begin();
                 id_itr != thread_data[i].dbNamesToHierachyIdentifierCounts[db_itr->first].end();
                 id_itr++) {
                 if ( writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first].find(id_itr->first) == writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first].end()) {
                     // create id slot if not present
                     writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first][id_itr->first] = 0;
                 }
                 // iterate global count
                 writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first][id_itr->first]++;
            }
        }
        
        // Writeout metacyc results to .pf file
        cout << "MetaCyc hits:" << thread_data[i].metaCycHits.size() << endl;
        
        for (vector<ANNOTATION>::iterator mc_itr = thread_data[i].metaCycHits.begin();
            mc_itr != thread_data[i].metaCycHits.end(); 
            mc_itr++) {
            
            end_base = start_base + mc_itr->length;
            
            sprintf(start_base_str, "%d", start_base);
            sprintf(end_base_str, "%d", end_base);
            writePfEntry(mc_itr->orf_id,
                         mc_itr->product,
                         mc_itr->ec,
                         string(start_base_str),
                         string(end_base_str), 
                         pf_output);
            // update startbase
            start_base = start_base + mc_itr->length + 10;
        }
        
        
        thread_data[i].clear();
        
    }
    pf_output.close();
    
    return (void *)NULL;

}

/*
 * Print out MetaCyc word trie using a stack-based recursive post-order traversal
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
    for (int i = 0; i < annotation_words.size(); i++) {
        word_list.clear();
        ptools_ptr = root; // reset to root
        for (int j = i; j < annotation_words.size(); j++) {
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
    pf_entry = pf_entry + "ID " + orf_id + "\n";
    pf_entry = pf_entry + "NAME " + orf_id + "\n";
    pf_entry = pf_entry + "STARTBASE " + start_base + "\n";
    pf_entry = pf_entry + "ENDBASE " + end_base + "\n";
    pf_entry = pf_entry + "PRODUCT-TYPE P" + "\n";
    pf_entry = pf_entry + "FUNCTION " + annotation_product + "\n"; // remove extract character
    if (ec_number != "")
        pf_entry = pf_entry + "EC " + ec_number + "\n";
    pf_entry = pf_entry + "//" + "\n";
    output << pf_entry;
}

void writePToolsAdminFiles(string ptools_dir, string sample_name) {

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
    genetic_elements_text = genetic_elements_text + "ID      0" + "\n";
    genetic_elements_text = genetic_elements_text + "NAME    0" + "\n";
    genetic_elements_text = genetic_elements_text + "TYPE    :READ/CONTIG" + "\n";
    genetic_elements_text = genetic_elements_text + "ANNOT-FILE      0.pf" + "\n";

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
    organism_params_text = organism_params_text + "ID      " + sample_name + "\n";
    organism_params_text = organism_params_text + "STORAGE FILE" + "\n";
    organism_params_text = organism_params_text + "NAME    " + sample_name + "\n";
    organism_params_text = organism_params_text + "ABBREV-NAME     " + sample_name + "\n";
    organism_params_text = organism_params_text + "STRAIN  1" + "\n";
    organism_params_text = organism_params_text + "RANK    |species|" + "\n";
    organism_params_text = organism_params_text + "NCBI-TAXON-ID   12908" + "\n";

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

