#include "process_blastoutput.h"





void *process_lines( void *_data) {
    THREAD_DATA *data = static_cast<THREAD_DATA *>(_data);
    vector<string>::iterator it;
    string orfid;
    char buf[1000];
    vector<char *> fields; 

    int work =0;
    for(it = data->lines.begin(); it != data->lines.end(); it++) {
       split(*it, fields, buf,'\t');
       if( fields.size()!=12) continue;
       orfid = ShortenORFId(fields[0]);
       work++;
    }
    std::cout << "workerd on " << work << std::endl;

    
    
    return (void *)NULL;

}


void create_threads(int num_threads, THREAD_DATA *thread_data) {
    pthread_t *threads;
    if((threads = (pthread_t *)malloc(sizeof(pthread_t)*num_threads))==0) {
      std::cout << "Error in allocating threads" << std::endl;
    }
    int rc;

     std::cout << " create threads \n";
    for(int i = 0; i < num_threads; i++) {
       if((rc = pthread_create(&threads[i], NULL, process_lines, (void *)(thread_data+i)))) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
       }
    }


    void *status;
    for(int i=0; i<num_threads; i++) {
      rc = pthread_join(threads[i], &status);
      if (rc) {
         printf("ERROR; return code from pthread_join() is %d\n", rc);
         exit(-1);
      }
        printf("Main: completed join with thread %d having a status   of %ld\n",i,(long)status);
     }

}


void process_blastoutput(const Options& options, const GLOBAL_PARAMS &params) {
    OutputParser parser(options, params);
    parser.initialize();

    THREAD_DATA *thread_data = new THREAD_DATA[options.num_threads];
    
    parser.initializeBatchReading();
    while(parser.readABatch()) {
       std::cout << " creating \n";
       parser.distributeInput(thread_data);

       std::cout << " create \n";

       create_threads(options.num_threads, thread_data);
   
    } 
    parser.closeBatchReading();


}
