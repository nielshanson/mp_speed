#include "process_blastoutput.h"




void process_blastoutput(const Options& options, const GLOBAL_PARAMS &params) {
    OutputParser parser(options, params);
    pthread_t *threads;
    if((threads = (pthread_t *)malloc(sizeof(pthread_t)*options.num_threads))==0) {
      std::cout << "Error in allocating threads" << std::endl;
    }
    parser.initialize();

    THREAD_DATA *thread_data = new THREAD_DATA[options.num_threads];
    
    parser.initializeBatchReading();
    while(parser.readABatch()) {
       parser.distributeInput(thread_data);
       //launch threads
       //write results
   
    } 
    parser.closeBatchReading();


}
