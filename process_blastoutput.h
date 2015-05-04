#ifndef _PARSE_BLASTOUTPUT
#define _PARSE_BLASTOUTPUT
#include <pthread.h>
#include <stdio.h>
#include <fstream>
#include <stdio.h>


#include "utilities.h"
#include "types.h"
#include "outputparser.h"
#include "externalsort.h"


void process_blastoutput(const Options &options, const GLOBAL_PARAMS &params);

#endif
