#include <map>
#include <vector>
#include <algorithm>
#include <assert.h>
#include "parse.h"

using namespace std;

GLOBAL_PARAMS params;

int main( int argc, char **argv ){
    // Parse options
    Options options;
    // if (argc < 9) { options.print_usage(argv[0]); exit(1); }
    if (options.SetOptions(argc, argv)==false) {
        // Something went wrong in arguments, print usage
        options.print_usage(argv[0]); exit(0);
    }
    if( !options.check_arguments() ) return 0;


    params.lambda = options.lambda; 
    params.k = options.k; 

    if( params.lambda==-1 or params.k ==-1) {
      if( options.algorithm==string("LAST")) {
          params.lambda = 0.300471;
          params.k = 0.103946;
      }
      if( options.algorithm==string("BLAST")) {
          params.lambda = 0.267;
          params.k = 0.0410;
      }
    }

    process_blastoutput(options, params);
    //options.print_options();

} 
