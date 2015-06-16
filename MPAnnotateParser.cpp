#include "MPAnnotateParser.h"


using namespace std;
#define PRINT_INTERVAL 5000000
#define BATCH_SIZE_PER_CORE 5000

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
    float evalue, prevevalue;
    
    if (this->options.debug) {
        std::cout << "input size " << i << "  " << this->inputbuffer.size() << std::endl;
    }
    
    prevorfid = "";
    prevevalue = 100;
    for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
        //  std::cout << db_info.db_names[i] << std::endl;
        int j = 0;
        // interate through the results from a given database and select the annotation with the
        // smallest evalue for each ORF
        prevevalue = 100;
        prevorfid = "";
        for(it =  dbwise_inputs[db_info.db_names[i]].begin(); it !=  dbwise_inputs[db_info.db_names[i]].end(); it++ ) {
            // for each annotation
            
            // extract ORF and evalue
            orfid = orf_extractor_from_blast((*it).c_str());
            evalue = evalue_extractor_from_blast((*it).c_str());

            // calculate thread to send result to
            bucketIndex = hashIntoBucket(orfid.c_str(), options.num_threads);
            
            // If database not in annot_objects, create it
            if(thread_data[bucketIndex].annot_objects.find(db_info.db_names[i])==thread_data[bucketIndex].annot_objects.end())
                thread_data[bucketIndex].annot_objects[db_info.db_names[i]]= map<string, ANNOTATION *>();
            
            // if this is new orfid and not first entry (header)
            if(prevorfid != orfid && j != 0 ) {
                prevevalue  = 100; // reset evalue
                // std::cout << orfid << "\t" << annotation->product << std::endl;
                thread_data[bucketIndex].orfids.push_back(orfid);
                annotation = createAnnotation((*it).c_str(), db_info.db_names[i]);
                thread_data[bucketIndex].annot_objects[db_info.db_names[i]][orfid] = annotation;
            }
            else if(prevevalue > evalue ) {
                annotation = createAnnotation((*it).c_str(), db_info.db_names[i]);
                thread_data[bucketIndex].annot_objects[db_info.db_names[i]][orfid] = annotation;
            }

            prevevalue = evalue; prevorfid = orfid;
            j++;
        }
    }

//  thread_data[bucketIndex].annot_objects.push_back(*it);
    if (this->options.debug) {
        for(unsigned int i=0; i< options.num_threads; i++) {
            std::cout << "Thread data " << i <<  std::endl;
            for(unsigned int j = 0; j < db_info.db_names.size(); j++ ) {
                std::cout << "  " << db_info.db_names[j] << " : " <<\
                thread_data[i].annot_objects[db_info.db_names[j]].size() << std::endl;
            }    
        }
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

    map<string, std::ifstream *>::iterator it;

    if(count == 0) return false;

    std::ifstream *parsedinput;
    for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
        if(leftover_lines[db_info.db_names[i]].size() > 0)
            dbwise_inputs[db_info.db_names[i]].push_back(leftover_lines[db_info.db_names[i]]);

        parsedinput = parsed_file_streams[db_info.db_names[i]];
        while( std::getline(*parsedinput, line ).good()) {
            if( line.size() > 0 && line[0]=='#') continue; // skip commented lines
//           std::cout << line << std::endl;
            orfid = orf_extractor_from_blast(line);
            if( gfforfs.find(orfid) == gfforfs.end()) {
                // save line
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

MPAnnotateParser::~MPAnnotateParser() {

}
