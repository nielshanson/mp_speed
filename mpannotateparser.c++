#include "MPAnnotateParser.h"


using namespace std;
#define PRINT_INTERVAL 5000000
#define BATCH_SIZE_PER_CORE 1000000


MPAnnotateParser::MPAnnotateParser(const MPAnnotateOptions &options,const DB_INFO & db_info){
    this->options = options;
    this->BATCH_SIZE = options.num_threads*BATCH_SIZE_PER_CORE;
    this->db_info= db_info;
}


void  MPAnnotateParser::closeBatchReading() {
   this->input.close();

   return ;
}


void  MPAnnotateParser::initializeBatchReading() {
   string filename  = options.input_gff;
   std::cout << " initialized batch reading gff file " << filename << std::endl;
   this->input.open(filename.c_str(), std::ifstream::in);
   if(!this->input.good()){
          std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
         return ;
   }

   std::ifstream *parsedinput ;
 
   for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
       std::ifstream *parsedinput  = new std::ifstream;;
       parsedinput->open((options.blast_dir + string("/") + db_info.input_blastouts[i]).c_str(), std::ifstream::in);
       parsed_file_streams[db_info.db_names[i]]  = parsedinput;
       leftover_lines[db_info.db_names[i]] = "";
       dbwise_inputs[db_info.db_names[i]] = vector<string>();
   } 


}

void MPAnnotateParser::distributeInput(THREAD_DATA_ANNOT *thread_data) {
    std::cout << "input size " << this->inputbuffer.size() << std::endl;
    string orfid;
    vector<char *>fields;
    int bucketIndex;

    for(unsigned int i=0; i< options.num_threads; i++) {
      thread_data[i].lines.clear();
    }
    char tempbuf[1000];


    ANNOTATION annotation;
/*
    for(vector<string>::iterator it = this->inputbuffer.begin(); it != this->inputbuffer.end(); it++) {
        orfid = orf_extractor_from_gff(*it);
        bucketIndex = hashIntoBucket(orfid.c_str(), options.num_threads); 
        split(*it, fields, buf,'\t');

        if( fields.size()!=9) {
            // Not a parsed blast file
            cerr << "Parsed BLAST/LASTout file  did not have the 10 columns." << endl;
            input.close();
            return ;
        }   

        annotation.bsr = atof(fields[4]);
        annotation.ec = fields[8]; 
        annotation.product = fields[9];

    }
*/

    std::cout << "input size " << this->inputbuffer.size() << std::endl;
    vector<string>::iterator it;

    for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
        std::cout << "input size " << i  << "  "  << this->inputbuffer.size() << std::endl;
        for(it =  dbwise_inputs[db_info.db_names[i]].begin(); it !=  dbwise_inputs[db_info.db_names[i]].end(); it++ ) {
           bucketIndex = hashIntoBucket((*it).c_str(), options.num_threads); 

           if( thread_data[bucketIndex].annot_objects.find(db_info.db_names[i])==thread_data[bucketIndex].annot_objects.end()) 
                thread_data[bucketIndex].annot_objects[db_info.db_names[i]]= map<string, ANNOTATION>();

           thread_data[bucketIndex].annot_objects[db_info.db_names[i]][*it] = annotation;
        }
        std::cout << " data for thread " <<   thread_data[bucketIndex].annot_objects[db_info.db_names[i]].size() << std::endl;
    }


     //   thread_data[bucketIndex].annot_objects.push_back(*it);
    for(unsigned int i=0; i< options.num_threads; i++) {
      std::cout << "Data size of thread " << i << "  " << thread_data[i].lines.size() << std::endl;
    }
    

}

bool MPAnnotateParser::readABatch() {
   std::cout << "Reading a new batch\n";
   int count = 0;
   string line;
   char buf[10000];
   vector<char *> fields;

   this->inputbuffer.clear();

  
   map<string, bool>  gfforfs;
   gfforfs.clear();

   for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
       dbwise_inputs[db_info.db_names[i]].clear();               
   }



   string orfid ;
   while( std::getline(this->input, line ).good()) {
        this->inputbuffer.push_back(line); 
        orfid = orf_extractor_from_gff(line); 
        gfforfs[orfid] = true;
        count++; 
        if(count > this->BATCH_SIZE - 1) { break;}
        if(count %PRINT_INTERVAL==0) 
          std::cout << count << std::endl;
   }
   
   std::cout << "orfs extracted " << gfforfs.size() << std::endl; 
   map<string, std::ifstream *>::iterator it;


   std::cout << "num databaessss  " << parsed_file_streams.size() << std::endl; 
   std::cout << "num databaessss  " << db_info.db_names.size() << std::endl; 

   if( count ==0) return false;


   std::ifstream *parsedinput;
   for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
       
       if(leftover_lines[db_info.db_names[i]].size() > 0)
           dbwise_inputs[db_info.db_names[i]].push_back(leftover_lines[db_info.db_names[i]]);               

       parsedinput = parsed_file_streams[db_info.db_names[i]];
       while( std::getline(*parsedinput, line ).good()) {
           if( line.size() > 0 && line[0]=='#') continue;
           orfid = orf_extractor_from_blast(line); 
           if( gfforfs.find(orfid) == gfforfs.end()) {
              leftover_lines[db_info.db_names[i]] = line;               
              break;
           }
           dbwise_inputs[db_info.db_names[i]].push_back(line);               
       }
   }

   std::cout << " ddatabaes one "<< db_info.db_names.size() << std::endl;;
   for(unsigned int i = 0; i < db_info.db_names.size(); i++ ) {
     std::cout << " db hits ->" << db_info.db_names[i] << "<=  "<< dbwise_inputs[db_info.db_names[i]].size() << std::endl;              
   }

   return true;;
   exit(0);

   



   if(count>0) return true; 

   return false;
}

MPAnnotateParser::~MPAnnotateParser() {

}
