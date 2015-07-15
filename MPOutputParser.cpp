#include "MPOutputParser.h"


using namespace std;
#define PRINT_INTERVAL 5000000
#define BATCH_SIZE_PER_CORE 1000000


OutputParser::OutputParser(const MPParseBlastOptions &options, const GLOBAL_PARAMS &params){
    this->options = options;
    this->params = params;
    this->ln2 =  0.69314718055994530941;
    this->lnk = log(params.k);
    this->BATCH_SIZE = options.num_threads*BATCH_SIZE_PER_CORE;

    this->r = new regex_t;
    this->regex_text = "([[:digit:]]+[_][[:digit:]]+)$";
    compile_regex(r, regex_text);

}



void OutputParser::initialize() {
    //string temp = options.input_blastout + ".tmp";
    // sort_parsed_blastouput(string("/tmp/"), options.input_blastout, temp,  1000);


//    std::cout << "Reading the query names " << std::endl;
    this->create_query_dictionary();
    //  std::cout <<  "Collected " << query_dictionary.size() << " unique  queries" << std::endl;
}



void OutputParser::create_query_dictionary() {
    string filename  = options.input_blastout;
    std::cout << "Filename : "  << filename <<  std::endl;
    this->input.open(filename.c_str(), std::ifstream::in);
    if(!this->input.good()){
        std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
        return ;
    }

    int count =0;
    string line;
    while( std::getline(this->input, line ).good()) {
        split(line, fields, this->buf,'\t');
        if( fields.size()!= 12) continue;
        query_dictionary[fields[1]] = true;
        if (count%PRINT_INTERVAL==0)
            std::cout << count << std::endl;
        count++;
    }

    this->input.close();
    std::cout <<  "Read " << count << " lines" << std::endl;
}




void OutputParser::create_annotation_dictionary( map<string, string> *annot_map ){

    string filename  = options.database_map;
    std::cout << "Reading annotation dictionary " <<  std::endl;
    std::cout << "File Name " <<  filename << std::endl;
    this->input.open(filename.c_str(), std::ifstream::in);
    if(!this->input.good()){
        std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
        return ;
    }

    int count =0;
    string line;
    string name, annotation;

    while( std::getline(this->input, line ).good()) {

        //    std::cout << line << std::endl;

        if( line.size()==0 or  line[0]!='>') continue;

        split_seq_name(line, fields, this->buf);

        name = (fields[0]);

        //  std::cout << name << std::endl;
        //std::cout << fields[1] << std::endl;

        if( fields.size()< 2)
            annotation = "hypothetical protein";
        else
            annotation = string(fields[1]);

        if( query_dictionary.find(name) != query_dictionary.end() )
            annot_map->insert(std::make_pair(name,annotation));

        if (count%PRINT_INTERVAL==0)
            std::cout << count << std::endl;
        count++;
    }

    std::cout << "Number of annotatons scanned " << count << std::endl;
    std::cout << "Number of annotation loaded " <<  annot_map->size() << std::endl;

    this->input.close();
}

void OutputParser::create_refBitScores(THREAD_DATA *thread_data) {

    string filename  = options.refscore_file;

    std::cout << "Reading refscores " <<  std::endl;
    std::cout << "Filename " << filename << "\n";
    this->input.open(filename.c_str(), std::ifstream::in);
    if(!this->input.good()){
        std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
        return ;
    }

    int count =0, bucketIndex;

    int *counts  = (int *)calloc(options.num_threads, sizeof(int));
    string line;
    string orfid;
    while( std::getline(this->input, line ).good()) {
        split(line, fields, this->buf,'\t');
        if( fields.size()!=2) continue;
        orfid = ShortenORFId(fields[0]) ;
        bucketIndex = hashIntoBucket(orfid.c_str(), options.num_threads);
        thread_data[bucketIndex].refscorePairs.push_back( std::make_pair<string, string>(orfid, fields[1]));
        thread_data[bucketIndex].refscores[orfid] = atof(fields[1]);
        // experiment add to all thread data instances
//        for (int i = 0; i < options.num_threads; i++) {
//            thread_data[i].refscorePairs.push_back( std::make_pair<string, string>(orfid, fields[1]));
//            thread_data[i].refscores[orfid] = atof(fields[1]);
//        }

        //refBitScores[orfid] = int((params.lambda*float(atof(fields[1])) - this->lnk )/this->ln2);

        counts[bucketIndex]++;
        if (count%PRINT_INTERVAL==0)
            std::cout << "x " << count << std::endl;
        count++;
    }
    this->input.close();

    for(unsigned int i =0; i < options.num_threads; i++)  {
        thread_data[i].ln2 = this->ln2;
        thread_data[i].lnk = this->lnk;
        thread_data[i].lambda = params.lambda;
    }

    std::cout << "Number of refscores loaded " <<  count << std::endl;


};


void  OutputParser::closeBatchReading() {
    this->input.close();

    return ;
}


void  OutputParser::initializeBatchReading() {
    string filename  = options.input_blastout;
    this->input.open(filename.c_str(), std::ifstream::in);
    if(!this->input.good()){
        std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
        return ;
    }

    return ;

}

void OutputParser::distributeInput(THREAD_DATA *thread_data) {
    string orfid;
    vector<string>::iterator it;
    int bucketIndex;


    for(unsigned int i=0; i< options.num_threads; i++) {
        thread_data[i].lines.clear();
    }

    std::cout << "input size " << this->inputbuffer.size() << std::endl;

    for(it = this->inputbuffer.begin(); it != this->inputbuffer.end(); it++) {

        split(*it, fields, this->buf,'\t');
        if( fields.size()!=12) continue;
        orfid = ShortenORFId(fields[0]) ;
        bucketIndex = hashIntoBucket(orfid.c_str(), options.num_threads);
        thread_data[bucketIndex].lines.push_back(*it);
    }

}
bool OutputParser::readABatch() {
    std::cout << "Reading a new batch\n";
    int count = 0;
    string line;

    this->inputbuffer.clear();

    while( std::getline(this->input, line ).good()) {
        this->inputbuffer.push_back(line);

        count++;
        if(count > this->BATCH_SIZE - 1) { return true;}

        if(count %PRINT_INTERVAL==0)
            std::cout << count << std::endl;
    }

    if(count>0) return true;

    return false;
}

OutputParser::~OutputParser() {

}
