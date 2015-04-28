#ifndef _PARSE_BLASTOUTPUT
#define _PARSE_BLASTOUTPUT
#include <pthread.h>
#include <stdio.h>

#include "utilities.h"
#include "types.h"
#include "outputparser.h"


void process_blastoutput(const Options &options, const GLOBAL_PARAMS &params);

#endif
