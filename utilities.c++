#include "utilities.h"


char *split_n_pick(const string  &strn,  char *buf, char d, unsigned int n) {
  strcpy(buf, strn.c_str());
   
  char *v=buf;
  char *s1 = buf;
  v=s1;

  unsigned int i =0; 

  while(*s1 != '\0') {
     if(*s1==d) { 
       *s1 = '\0';
       i++;
       if(i > n) return  v ;
       v=s1+1;
     }
     s1++;
  }
  return v;
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

void split_seq_name(const string  &strn, std::vector<char *> &v, char *buf) {
    strcpy(buf, strn.c_str());

    if(buf[0]!='>') {
        v.push_back(buf);
        return;
    }

    char *s1 = buf+1;
    v.clear();
    v.push_back(s1);

    while(*s1 != '\0') {
        if(*s1==' ') {
            *s1 = '\0';
            v.push_back(s1+1);
            break;
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

    if( !fromstart && pos > 0) return true;
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
            p++;
            //*p++;
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

string ShortenORFId(const string &s) {

    const char * p = s.c_str();
    char BUFFER[200];
    strcpy(BUFFER, p);

    char *c;

    c = BUFFER ;

    c = c + strlen(BUFFER)-1;

    unsigned int S =0;

    while( c >= BUFFER ) {
        if(S==0) {
            if( isdigit(*c) || *c=='_') {
                if( *c=='_') S++;
            }
            else{
                return string(BUFFER);
            }
        }
        else if(S==1) {
            if( isdigit(*c))
                S++;
            else
                return string(BUFFER);
        }
        else {
            if( !isdigit(*c) ) {
                break;
            }
        }
        c--;
    }

    return string(c+1);
}


int compile_regex(regex_t * r, const char * regex_text) {
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

string getpattern(regex_t *r , const char *to_match, unsigned int no ) {
    /* "P" is a pointer into the string which points to the end of the previous match. */
    const char * p = to_match;
    /* "N_matches" is the maximum number of matches allowed. */ /* "M" contains the matches found. */
    regmatch_t m[100];
    char buf[1000];
    int i = no;
    int nomatch = regexec(r, p, no+1, m, 0);
    if (nomatch)  return string();
    int start;
    int finish;
    if (m[i].rm_so == -1) return string();
    start = m[i].rm_so + (p - to_match);
    finish = m[i].rm_eo + (p - to_match);
    memcpy(buf, to_match+ start, finish-start);
    buf[finish-start]='\0';
    return string(buf);
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

/*
 * Process the product field of the BLAST/LASTout.parsed file to extract taxonomy
 * contained in square brackets '[' ']'. Returns the expected taxonomy if found. Returns
 * 'no-taxonomy' otherwise.
 */
string getTaxonomyFromProduct(const char *str) {
    char buf[10000];
    const char *c, *b, *e;
    bool found = false;
    c = str;

    bool front = false;
// TODO: Basic implementaiton. Occassionally will run into problems with double taxonomies and '[[' ']]'
    while( *c!='\0') {
        // continue iterating through string until end
        if(*c=='[') {
            front = true;
            b = c+1;
            c++;
            continue;
        }
        if(*c==']' && front) {
            e = c;
            found = true;
        }
        c++;
    }

    if( found ) {
        unsigned int len = e - b;
        memcpy(buf,b, len);
        buf[len] ='\0';
        return string(buf);
    }

    return "no-taxonomy";
}

string getECNo(const char *str, unsigned int d) {
    char buf[100];
    unsigned int S=0;
    const char *c, *b, *e;
    bool found = false;

    c = str;
    while( *c!='\0') {

        if(S%2==0)  {
            if( isdigit(*c)) {
                if(S==0) b = c;
                S++;
                if(S==2*d+1) {
                    found = true;
                    e = c+1;
                }
            }
            else
                S=0;
        }
        else {// S%2 ==1
            if(S < 2*d +1 ) {
                if( isdigit(*c))
                    ;
                else if(*c=='.')
                    S++;
                else
                    S=0;
            }
            else  {//(S = 2*d +1 )
                if(!isdigit(*c)) {
                    e = c;
                    break;
                }
                else{
                    e = c+1;
                }

            }
        }
        c++;
    }


    if( found ) {
        unsigned int len = e - b;
        memcpy(buf,b, len);
        buf[len] ='\0';
        return string(buf);
    }

    buf[0]='\0';
    return string(buf);
}
