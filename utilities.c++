#include "utilities.h"
//#include <iostream>
//#include <fstream>


void Options::print_usage(char *arg) {
   std::cout << "USAGE : "   << arg <<"\n"\
             << "      : -b  <blastoutput>    [REQUIRED]\n"\
             << "      : -d  <database name>  [REQUIRED]\n"\
             << "      : -o  <parsedoutput>   [REQUIRED]\n"\
             << "      : -r  <refscorefile>   [REQUIRED]\n"
             << "      : -m  <databaemap>     [REQUIRED]\n"\
             << "      : -a  <algorihm>       [OPTIONAL,     default: BLAST]\n"\
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

 

void split(const string  &strn, std::vector<char *> &v, char *buf, char d) {
  strcpy(buf, strn.c_str());
  char *s1 = buf;
  v.clear();
  v.push_back(s1);
  while(*s1 != '\0') {
     if(*s1==d) { 
       *s1 = '\0';
       v.push_back(s1+1);
     }
     s1++;
  }
}

std::string get_orf_name(std::string  &strn, std::vector<char *> &v, char *buf) {
    split(strn, v, buf, ';'); 
    if(v.size() == 0)  return std::string("");
    split(std::string(v[0]), v, buf, '=');
    if(v.size() < 2)  return std::string("");
    return std::string(v[1]);
}

bool matchString(const string &str, const string & stringtomatch, bool fromstart) {

    unsigned long pos = str.find(stringtomatch);
    if(fromstart &&  pos ==0 ) return true;

    if( !fromstart && pos >= 0) return true;
    return false;

}

string to_string(unsigned long i) {
    char  c[100];
    char *p = c;
    int j = 0;
    char z = '0';

    while( i > 0 ) {
       if(i< 10) {
         *p='0' + i;
          p++;
          break;
       } 
       else {
           j = i%10;
           i = (i - j)/10;
           *p =  z + j;
           *p++;
        }
    }
    *p = '\0';
    p--;

    return string(c);
}


char BUFFER[1000];
string ShortenORFId(const string &s, regex_t *r) {
    const char * p = s.c_str();
    char *buf = BUFFER;
    regmatch_t m[1];
    int nomatch = regexec(r, p, 1, m, 0);
    if (nomatch) {
        return s;
    }
    p += m[0].rm_so;
    int d = m[0].rm_eo - m[0].rm_so;
    while( d > 0) {
      *buf = *p;
      d--;
      buf++;  p++;
    }
    *buf = '\0';
    return string(BUFFER);
}


int compile_regex(regex_t * r, const char * regex_text)
{
    int status = regcomp(r, regex_text, REG_EXTENDED|REG_NEWLINE);
    if (status != 0) {
    char error_message[MAX_ERROR_MSG];
    regerror (status, r, error_message, MAX_ERROR_MSG);
        printf ("Regex error compiling '%s': %s\n",
                 regex_text, error_message);
        return 1;
    }
    return 0;
}


int hashIntoBucket(const char *str, unsigned int index) {
    int hashValue = 0;
    char buffer[1000];

    char *first, *second;
    first  = buffer;

    char *x = buffer;
    const char *p =str; 
    
    // Extract contig id and orf_id from string
    while( *p != '\0') {
      if( *p=='_' ) {
         *x = '\0';
         second = x +1;
      }
      else 
         *x = *p;

      p++; 
      x++;
    }
    *x = '\0';

    //std::cout << str <<  "  " << first << "  " << second << std::endl;
    
    // Find the longer of the two ids
    char *dest, *destfixed,  *src; 
    dest = first;
    src = second;
    int lenextra = strlen(second) - strlen(first);
    if( lenextra > 0)  {
      dest = second;
      src = first;
    }
    else
        lenextra = -lenextra;

    destfixed = dest; // always remember the initial point of longer

    dest = dest + lenextra; // move to aligned position
    
    // XOR aligned bits of source and destination starting at aligned position
    while(*src != '\0') {
       *dest = (*dest ) ^ (*src);
       dest++; src++;
    }
    
    // feed into uniform hasing algorithm
    while( *destfixed != '\0') {
       //hashValue = (hashValue << 4) + (unsigned int)(*destfixed);
       hashValue = hashValue + (unsigned int)(*destfixed);
/*
       int hiBits = hashValue  & 0xF0000000;
       if(hiBits!=0) 
          hashValue ^= hiBits>> 24; 
       hashValue &= ~hiBits;
*/
       destfixed++;
    }   
    return hashValue%index;
}

