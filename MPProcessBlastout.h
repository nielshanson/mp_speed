#ifndef _PARSE_BLASTOUTPUT
#define _PARSE_BLASTOUTPUT
#include <pthread.h>
#include <stdio.h>
#include <fstream>
#include <stdio.h>

#include "utilities.h"
#include "types.h"
#include "MPOutputParser.h"
#include "externalsort.h"

void MPProcessBlastout(const MPParseBlastOptions &options, const GLOBAL_PARAMS &params);

#endif
