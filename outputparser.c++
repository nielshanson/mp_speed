#include "outputparser.h"


using namespace std;



OutputParser::OutputParser(const Options &options, const GLOBAL_PARAMS &params){
    this->options = options;
    this->params = params;
    this->ln2 =  0.69314718055994530941;
    this->lnk = log(params.k);
    this->BATCH_SIZE = options.num_threads*1000000;

    this->r = new regex_t;
    this->regex_text = "([[:digit:]]+[_][[:digit:]]+)$";
    compile_regex(r, regex_text);

}






void OutputParser::initialize() {
 //   this->create_query_dictionary();
    std::cout << "done  create query " << std::endl;
//    this->create_annotation_dictionary();
    std::cout << "done  create dciotnary query " << std::endl;
//    this->create_refBitScores() ;
    std::cout << "done  refscore " << std::endl;

}



void OutputParser::create_query_dictionary() {
   string filename  = options.input_blastout;
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
      if (count%1000000==0) 
       std::cout << count << std::endl;
      count++;
   }   
       
   this->input.close();
}




void OutputParser::create_annotation_dictionary( map<string, string> &annot_map ){

   string filename  = options.database_map;
   this->input.open(filename.c_str(), std::ifstream::in);
   if(!this->input.good()){
          std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
         return ;
   }

   int count =0;
   string line;
   string name, annotation;

   while( std::getline(this->input, line ).good()) {

      if( line.size()==0 or  line[0]!='>') continue;

      split(line, fields, this->buf,'\t');
      name = (fields[0] +1);

      if( fields.size()< 2) 
        annotation = "hypothetical protein";
      else
        annotation = string(fields[1]);

      annot_map[name] = annotation;
      if (count%1000000==0) 
         std::cout << count << std::endl;
      count++;
   }   
       
   this->input.close();




}
void OutputParser::create_refBitScores(THREAD_DATA *thread_data) {

   string filename  = options.refscore_file;
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
      //refBitScores[orfid] = int((params.lambda*float(atof(fields[1])) - this->lnk )/this->ln2);

      counts[bucketIndex]++;
      if (count%1000000==0) 
         std::cout << "x " << count << std::endl;
      count++;
   }   
   this->input.close();


   for(unsigned int i =0; i < options.num_threads; i++)  {
       thread_data->ln2 = this->ln2;
       thread_data->lnk = this->lnk;
       thread_data->lambda = params.lambda;
   }

/*
   std::cout << "done " << std::endl;
   for( int i =0; i < options.num_threads; i++) 
      std::cout << i << "  "<< counts[i] << std::endl;
*/
    


       
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

    int *counts  = (int *)calloc(options.num_threads, sizeof(int));

    std::cout << "input size " << this->inputbuffer.size() << std::endl;

    for(it = this->inputbuffer.begin(); it != this->inputbuffer.end(); it++) {
         
      split(*it, fields, this->buf,'\t');
      if( fields.size()!=12) continue;
      orfid = ShortenORFId(fields[0]) ;
      bucketIndex = hashIntoBucket(orfid.c_str(), options.num_threads); 
      thread_data[bucketIndex].lines.push_back(*it);
  //    counts[bucketIndex]++;
    }



/*
    std::cout << "done " << std::endl;
    for( int i =0; i < options.num_threads; i++) 
       std::cout << i << "  "<< counts[i] << std::endl;
*/


}
bool OutputParser::readABatch() {
   std::cout << "Reading a new batch\n";
   int bucketIndex ; 
   int count = 0;
   string line;

   this->inputbuffer.clear();

   while( std::getline(this->input, line ).good()) {
      this->inputbuffer.push_back(line); 

      count++;
      if(count > this->BATCH_SIZE - 1) { return true;}

      if(count %1000000==0) 
        std::cout << count << std::endl;
   }

   if(count>0) return true; 

   return false;
}

OutputParser::~OutputParser() {

}
