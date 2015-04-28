#include "outputparser.h"


using namespace std;



OutputParser::OutputParser(const Options &options, const GLOBAL_PARAMS &params){
    this->options = options;
    this->params = params;
    this->ln2 =  0.69314718055994530941;
    this->lnk = log(params.k);

}


void OutputParser::initialize() {
    this->create_query_dictionary();
    std::cout << "done  create query " << std::endl;
    this->create_annotation_dictionary();
    std::cout << "done  create dciotnary query " << std::endl;
    this->create_refBitScores() ;
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




void OutputParser::create_annotation_dictionary(){

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
void OutputParser::create_refBitScores() {

   string filename  = options.refscore_file;
   this->input.open(filename.c_str(), std::ifstream::in);
   if(!this->input.good()){
          std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
         return ;
   }

   int count =0;
   string line;
   while( std::getline(this->input, line ).good()) {
      split(line, fields, this->buf,'\t');
      if( fields.size()!=2) continue;
      refBitScores[fields[0]] = int((params.lambda*float(atof(fields[1])) - this->lnk )/this->ln2);
      if (count%1000000==0) 
         std::cout << count << std::endl;
      count++;
   }   
       


};


OutputParser::~OutputParser() {


}
