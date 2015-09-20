#include "MPAnnotateParser.h"


using namespace std;
#define PRINT_INTERVAL 5000000
#define BATCH_SIZE_PER_CORE 100000

/*
 * Initialize MPAnnotateParser. Set batch size to the number of threads * core batch size
 */
MPAnnotateParser::MPAnnotateParser(const MPAnnotateOptions &options, const DB_INFO & db_info){
    this->options = options;
    this->BATCH_SIZE = options.num_threads*BATCH_SIZE_PER_CORE;
    this->db_info = db_info;
}

void  MPAnnotateParser::closeBatchReading() {
    this->input.close();
    return ;
}

/*
 * Initializes the data structures for parallel parsing of BLASTout.parsed.txt. Opens a file pointer to the gff file
 * (sorted) for each thread. Creates an array of file-pointers, an array of left-over lines for spillover, and an
 * array of string vectors that collect the parsed data.
 */
void MPAnnotateParser::initializeBatchReading() {
    string filename  = options.input_gff;
    if (options.debug) {
        std::cout << "Initialized batch reading gff file " << filename << std::endl;
    }
    this->input.open(filename.c_str(), std::ifstream::in);
    if(!this->input.good()){
        std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
        return ;
    }
    
    std::ifstream *parsedinput ;
    for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
        // Open B/LASTout file
        parsedinput = new std::ifstream;
        parsedinput->open((options.blast_dir + string("/") + db_info.input_blastouts[i]).c_str(), std::ifstream::in);
        // Save fp into parsed_file_streams
        parsed_file_streams[db_info.db_names[i]]  = parsedinput;
        // Extra object for remaining lines if not within BATCH chunk
        leftover_lines[db_info.db_names[i]] = "";
        dbwise_inputs[db_info.db_names[i]] = vector<string>();
    }
}

/*
 * Split input data to be handled by threads
 */
void MPAnnotateParser::distributeInput(THREAD_DATA_ANNOT *thread_data) {
    vector<char *>fields;
    int bucketIndex;
    
    // clear out old data
    unsigned int i;
    for(i=0; i< options.num_threads; i++) {
        thread_data[i].lines.clear();
    }
    
    ANNOTATION *annotation;

    vector<string>::iterator it;
    string orfid, prevorfid;
    // float evalue, prevevalue;
    
    if (this->options.debug) {
        std::cout << "distributeInput()" << endl;
        std::cout << "Threads: " << i << " Parsed_GFF_Lines:" << this->inputbuffer.size() << std::endl;
    }

    // clear old annotations
    for (unsigned int k = 0; k < options.num_threads; k++) {
            thread_data[k].annot_objects.clear();
            thread_data[k].orfids.clear();
    }

    prevorfid = "";
    //prevevalue = 100;
    int max_num_hits_db = 5; // top hits pull from options
    int num_hits_db = 0;
    for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
        //  std::cout << db_info.db_names[i] << std::endl;
        int j = 0;
        // interate through the results from a given database and select the annotation with the
        // smallest evalue for each ORF
        // prevevalue = 100;
        prevorfid = "";
        num_hits_db = 0;
        for(it =  dbwise_inputs[db_info.db_names[i]].begin(); it !=  dbwise_inputs[db_info.db_names[i]].end(); it++ ) {
            // for each annotation
            // extract ORF and evalue
            orfid = orf_extractor_from_blast((*it).c_str());
            // evalue = evalue_extractor_from_blast((*it).c_str());

            // calculate thread to send result to
            bucketIndex = hashIntoBucket(orfid.c_str(), options.num_threads);
            
            // If database not in annot_objects, create it
            if(thread_data[bucketIndex].annot_objects.find(db_info.db_names[i])==thread_data[bucketIndex].annot_objects.end())
                thread_data[bucketIndex].annot_objects[db_info.db_names[i]]= map<string, vector< ANNOTATION *> >();

            if(prevorfid != orfid) {
                num_hits_db = 0;
                // std::cout << orfid << "\t" << annotation->product << std::endl;
                thread_data[bucketIndex].orfids.push_back(orfid);
            }
            if (num_hits_db < max_num_hits_db) {
                annotation = createAnnotation((*it).c_str(), db_info.db_names[i]);
                thread_data[bucketIndex].annot_objects[db_info.db_names[i]][orfid].push_back(annotation);
            }
            prevorfid = orfid;
            num_hits_db++;
            j++;
        }
        cout << db_info.db_names[i] << ":" << j << endl;
    }

//  thread_data[bucketIndex].annot_objects.push_back(*it);
    int grand_total_orfs = 0;
    int grand_total_anno = 0;
    int total_orfs = 0;
    int total_num_anno = 0;
    int orfs = 0;
    int num_anno = 0;
    if (this->options.debug) {
        for(unsigned int i=0; i< options.num_threads; i++) {
            std::cout << "Thread data " << i <<  std::endl;
            total_orfs = 0;
            total_num_anno = 0;
            for(unsigned int j = 0; j < db_info.db_names.size(); j++ ) {
                orfs = 0;
                num_anno = 0;
                orfs = thread_data[i].annot_objects[db_info.db_names[j]].size();
                for (map<string, vector< ANNOTATION *> >::iterator itr = thread_data[i].annot_objects[db_info.db_names[j]].begin();
                     itr != thread_data[i].annot_objects[db_info.db_names[j]].end();
                     ++itr) {
                     num_anno += (itr->second).size();
                }
                std::cout << "  " << db_info.db_names[j] << " : " << orfs << " ORFs with " << num_anno << " annotations" << std::endl;
                total_orfs += orfs;
                total_num_anno += num_anno;
            }
            std::cout << " Total ORFs " << total_orfs << " with " << total_num_anno << " annotations" << endl;
            grand_total_orfs += total_orfs;
            grand_total_anno += total_num_anno;
        }
        std::cout << "Grand Total ORFs " << grand_total_orfs << " with " << grand_total_anno << " annotations" << endl;
    }
    if (this-options.debug)
        cout << "End of distributeInput()" << endl;
}

/*
 * Reads a batch
 */
bool MPAnnotateParser::readBatch() {
    int count = 0;
    string line;
    vector<char *> fields;

    this->inputbuffer.clear();

    map<string, bool>  gfforfs;
    gfforfs.clear();

    for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
        dbwise_inputs[db_info.db_names[i]].clear();
    }

    // Read ORF_IDs from gff file to handle
    string orfid;
    while( std::getline(this->input, line ).good()) {
        this->inputbuffer.push_back(line);
        orfid = orf_extractor_from_gff(line);
        gfforfs[orfid] = true;
        count++;
        if(count > this->BATCH_SIZE - 1) { break; }
        if(count %PRINT_INTERVAL==0)
            std::cout << count << std::endl;
    }

    if(count == 0) return false;


    map<string, std::ifstream *>::iterator it;
    std::ifstream *parsedinput;
    for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
        if(leftover_lines[db_info.db_names[i]].size() > 0)
            dbwise_inputs[db_info.db_names[i]].push_back(leftover_lines[db_info.db_names[i]]);
            leftover_lines[db_info.db_names[i]].clear();

        parsedinput = parsed_file_streams[db_info.db_names[i]];
        while( std::getline(*parsedinput, line ).good()) {
            if( line.size() > 0 && line[0]=='#') continue; // skip commented lines
//           std::cout << line << std::endl;
            orfid = orf_extractor_from_blast(line);
            if( gfforfs.find(orfid) == gfforfs.end()) {
                // found ORF that was not in gff BATCH
                leftover_lines[db_info.db_names[i]] = line;
                break;
            }
            // add line
            dbwise_inputs[db_info.db_names[i]].push_back(line);
        }
    }

   if (options.debug) {
       cout << "readBatch(): " << endl;
       for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
           std::cout << "\t" << db_info.db_names[i] << ": " << dbwise_inputs[db_info.db_names[i]].size() << std::endl; 
       }
   }

    if(count>0) return true;

    return false;
}


/*
 * Takes a ncbi taxonomy id and returns it full lineage separated by semicolons
 */
string MPAnnotateParser::prepareRefSeqTaxonomy(string ncbi_id, map<string, string> NCBI_ID_to_Common, NCBITree* ncbi_tree) {
    string result = "";
    vector<string> id_lineage = ncbi_tree->getLineage(ncbi_id);
    reverse(id_lineage.begin(),id_lineage.end());
    result = NCBI_ID_to_Common[ncbi_id];
    for (unsigned int i = 1; i < id_lineage.size(); ++i) {
        result = result + ";" + NCBI_ID_to_Common[id_lineage[i]];
    }
    return result;
}



/*
 * Writes functional annotation hits to results folder as .tree.count.txt files
 */
void MPAnnotateParser::writeFunctionalHierarchyFiles(WRITER_DATA_ANNOT *writer_data, MPAnnotateOptions options) {
    string output_dir = options.results_dir;
    string db_name = "";
    string ending = ".tree.count.txt";
    string filename = "";
    string header = "ID\talias\tcount";
    string sample_name = options.sample_name;
    string alias_name = "";
    map<string, string> db_hierarchy_id_map;
    
    int total = 0;
    
    if (options.debug) cout << "globalDbNamesToHierachyIdentifierCounts:" << endl;
    
    for ( map<string, map<string, int> >::iterator db_itr = writer_data->globalDbNamesToHierachyIdentifierCounts.begin();
        db_itr != writer_data->globalDbNamesToHierachyIdentifierCounts.end();
        db_itr ++ ) {
            
        db_name = db_itr->first;
        string temp_db_str = "";
        temp_db_str = to_upper(db_name);

        if (writer_data->dbNamesToHierarchyIdentifierMaps.find(db_name) != writer_data->dbNamesToHierarchyIdentifierMaps.end()) {
            db_hierarchy_id_map = writer_data->dbNamesToHierarchyIdentifierMaps[db_name];
        } else if ( temp_db_str.find("REFSEQ") != std::string::npos ) {
            db_hierarchy_id_map = writer_data->ncbi_tree->NCBI_ID_to_Common;
        } else {
            db_hierarchy_id_map = map<string, string>();
        }
        
        if (writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first].size() >= 0) {
            
            filename = output_dir + "/" + sample_name + "." + db_name + ending;
            this->output.open(filename.c_str(), std::ofstream::out);
            
            if(!this->output.good()){
                std::cerr << "Error opening '" << filename <<"'. Bailing out." << std::endl;
                exit(-1);
            }
            
            // Write header
            output << header << endl;
            
            total = 0;
            // Iterate through ids
            for (map<string, int>::iterator id_itr = writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first].begin();
                      id_itr != writer_data->globalDbNamesToHierachyIdentifierCounts[db_itr->first].end();
                      id_itr++) {
                total += id_itr->second;
                if (db_hierarchy_id_map.size() > 0) {
                    if (db_hierarchy_id_map.find(id_itr->first) != db_hierarchy_id_map.end()) {
                        if (temp_db_str.find("REFSEQ") != std::string::npos) {
                            alias_name = prepareRefSeqTaxonomy(id_itr->first, db_hierarchy_id_map, writer_data->ncbi_tree);
                        } else {
                            alias_name = db_hierarchy_id_map[id_itr->first];
                        }
                    }
                }
                output << id_itr->first << "\t" << alias_name << "\t" << id_itr->second << endl;
                
            }
            output.close();
            
            if (options.debug) cout << sample_name + "." + db_name + ending << ": " << total << " functional hierachy IDs" << endl;
        } else {
            if (options.debug) cout << "Warning: " << sample_name + "." + db_name + ending << " contains no functional hierachy IDs!" << endl;
        }
    }
}

MPAnnotateParser::~MPAnnotateParser() {

}
