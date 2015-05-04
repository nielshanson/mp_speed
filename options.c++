#include "options.h"
//#include <iostream>
//#include <fstream>


void Options::print_usage(char *arg) {
   std::cout << "USAGE : "   << arg <<"\n"\
             << "      : -b  <blastoutput>    [REQUIRED]\n"\
             << "      : -d  <database name>  [REQUIRED]\n"\
             << "      : -o  <parsedoutput>   [REQUIRED]\n"\
             << "      : -r  <refscorefile>   [REQUIRED]\n"
             << "      : -m  <databaemap>     [REQUIRED]\n"\
             << "      : -a/--algorithm  <algorihm>       [OPTIONAL,     default: BLAST]\n"\
             << "      : -s  <sample_name>    [OPTIONAL,     default: BLAST]\n"\
             << "      : --min_score          <min_score>,   default: 20 \n"\
             << "      : --max_evalue         <maxEvalue>    [OPTIONAL, default 1e-6]\n"\
             << "      : --min_bsr            <minBSR>       [OPTIONAL, default 0.4]\n"\
             << "      : --min_length         <minlength>    [OPTIONAL, default 30]\n"\
             << "      : --min_identity       <min_identity> [OPTIONAL, default 30]\n"\
             << "      : --limit              <limit>        [OPTIONAL, deafult 5]\n"\
             << "      : --lambda             <lambda>       [OPTIONAL]\n"\
             << "      : --k                  <k>            [OPTIONAL]\n"\
             << "      : --num_threads        <t>            [OPTIONAL, default 1]\n"\
             << std::endl;
}

bool Options::SetOptions(int argc, char *argv[]) { 
   for(int i = 1; i < argc ; i++) {   
       if( strncmp(argv[i], "-b", strlen("-b")) == 0 ) {   
          this->input_blastout = argv[++i];
       }   
       else if( strncmp(argv[i], "-d", strlen("-d")) == 0 ) {   
          this->database_name = argv[++i];
       }   
       else if( strncmp(argv[i], "-o", strlen("-o")) == 0 ) {   
          this->parsed_output = argv[++i];
       }   
       else if(strncmp(argv[i], "-r", strlen("-r")) == 0 ) {   
          this->refscore_file =argv[++i];
       }   
       else if(strncmp(argv[i], "-m", strlen("-m")) == 0 ) {   
          this->database_map =argv[++i];
       }   
       else if( strncmp(argv[i], "-a", strlen("-a")) == 0 ) {   
          this->algorithm = argv[++i];
       }
       else if( strncmp(argv[i], "--algorithm", strlen("--algorithm")) == 0 ) {   
          this->algorithm = argv[++i];
       }
       else if( strncmp(argv[i], "-s", strlen("-s")) == 0 ) {   
          this->sample_name = argv[++i];
       }
       else if( strncmp(argv[i], "--min_score", strlen("--min_score")) == 0 ) {   
          this->min_score = atof(argv[++i]);
       }
       else if( strncmp(argv[i], "--min_bsr", strlen("--min_bsr")) == 0 ) {   
          this->min_bsr = atof(argv[++i]);
       }
       else if( strncmp(argv[i], "--max_evalue", strlen("--max_evalue")) == 0 ) {   
          this->max_evalue = atof(argv[++i]);
       }
       else if( strncmp(argv[i], "--min_length", strlen("--min_length")) == 0 ) {   
          this->min_length = atoi(argv[++i]);
       }
       else if( strncmp(argv[i], "--num_threads", strlen("--num_threads")) == 0 ) {   
          this->num_threads = atoi(argv[++i]);
       }
       else if( strncmp(argv[i], "--min_identity", strlen("--min_identity")) == 0 ) {   
          this->min_identity = atof( argv[++i]);
       }
       else if( strncmp(argv[i], "--limit", strlen("--limit")) == 0  ) {   
          this->limit = atoi(argv[++i]);
       }
       else if( strncmp(argv[i], "--lambda", strlen("--lambda")) == 0  ) {   
          this->lambda = atof(argv[++i]);
       }
       else if( strncmp(argv[i], "--k", strlen("--k")) == 0  ) {   
          this->k = atof(argv[++i]);
       }
       else if( strncmp(argv[i], "--k", strlen("--k")) == 0  ) {   
          this->k = atof(argv[++i]);
       }
       else {
            cout << "ERROR: Cannot recognize argument " << argv[i] << std::endl;;
            return false;
       }
    } //for loop for arguments processing

    return true;
    
};


bool Options::check_arguments(){

    if(this->input_blastout.size() == 0) {
         cout << "There sould be at least one blastoutput file"  << std::endl; 
         return false;
    }

    if(this->database_name.size() == 0) {
         cout << "There sould be at least one database name"   << std::endl;
         return false;
    }

    if(this->database_map.size()== 0) {
         cout << "There sould be at least one database map file name"  << std::endl;
         return false;
    }

    if(this->refscore_file.size()==0) {
       cout << "Must specify the refscore" << std::endl;
       return false;
    }

    return true;
}



void Options::print_options() {
}

