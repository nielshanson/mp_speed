#include "process_blastoutput.h"



void *write_results(void *writer_data) ;


void *compute_refscores( void *_data) {
    THREAD_DATA *data = static_cast<THREAD_DATA *>(_data);
    vector< std::pair<string, string> >::iterator it;
    string orfid;
    vector<char *> fields; 

    for(it = data->refscorePairs.begin(); it != data->refscorePairs.end(); it++) {
       data->refscores[it->first] =(data->lambda*atof(it->second.c_str())-data->lnk)/data->ln2;
    }

    data->refscorePairs.clear();
    return (void *)NULL;

}


bool isWithinCutoffs(const string &orfid, const vector<char *> &fields, BLASTOUT_DATA &data, THREAD_DATA *thread_data) {

  
  try{
     data.query = orfid; 
     data.target = fields[1];
     data.q_length = atoi(fields[7]) - atoi(fields[6]) + 1;
     data.bitscore = atof(fields[11]);
     data.bsr = atof(fields[11])/thread_data->refscores[orfid];
     data.expect = atof(fields[10]);
     data.aln_length = atof(fields[3]);
     data.identity = atof(fields[2]);
     if( thread_data->annot_map->find(fields[1])!=thread_data->annot_map->end())
        data.product=thread_data->annot_map->find(fields[1])->second;
     else
        data.product= "hypothetical protein";

     data.ec = getECNo(data.product.c_str(), 3);
     //std::cout << data.product << std::endl;
    // std::cout << "<<" << data.ec  << ">>"  << std::endl;


//    std::cout << "passed 0" << std::endl;  
    if(data.q_length < thread_data->options.min_length) return false;
    if(data.bitscore < thread_data->options.min_score) return false;
 //   std::cout << "passed 1" << std::endl;  
  //  std::cout <<  data.expect << "   " <<  thread_data->options.max_evalue << std::endl;
    if(data.expect > thread_data->options.max_evalue) return false;
   // std::cout << "passed 2" << std::endl;  
    if(data.identity < thread_data->options.min_identity) return false;
    if(data.bsr < thread_data->options.min_bsr) return false;

   }

   catch(...) {
      return false;
   }
   
   return true;
}



void *process_lines( void *_data) {
    THREAD_DATA *data = static_cast<THREAD_DATA *>(_data);
    vector<string>::iterator it;
    string orfid;
    char buf[1000];
    vector<char *> fields; 
    map<string, unsigned int> hit_counts;
    BLASTOUT_DATA blastout_data;

    int work =0;
    for(it = data->lines.begin(); it != data->lines.end(); it++) {
       split(*it, fields, buf,'\t');
       if( fields.size()!=12) continue;
       orfid = ShortenORFId(fields[0]);

       if(hit_counts.find(orfid) == hit_counts.end() )
          hit_counts[orfid] = 0;

       if(hit_counts[orfid] >= data->options.limit) 
           continue;

       blastout_data.setDefault();

       if( !isWithinCutoffs(orfid,  fields, blastout_data, data ))
         continue;

       data->output[data->b].push_back(blastout_data);

       hit_counts[orfid]++;

       work++;
    }

    data->lines.clear();


    //std::cout << "worked on " << work << std::endl;

    return (void *)NULL;

}

void create_threads_refscores(int num_threads, THREAD_DATA *thread_data ) {
    pthread_t *threads;
    if((threads = (pthread_t *)malloc(sizeof(pthread_t)*num_threads))==0) {
      std::cout << "Error in allocating threads" << std::endl;
    }
    int rc;

    for(int i = 0; i < num_threads; i++) {
       if((rc = pthread_create(&threads[i], NULL, compute_refscores, (void *)(thread_data+i)))) {
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
     }

}



void create_threads_parse(int num_threads, THREAD_DATA *thread_data, WRITER_DATA *writer_data) {
    pthread_t *threads;
    if((threads = (pthread_t *)malloc(sizeof(pthread_t)*num_threads))==0) {
      std::cout << "Error in allocating threads" << std::endl;
    }
    int rc;

    std::cout << "Number of threads created to process lines " << num_threads << std::endl;

    for(int i = 0; i < num_threads; i++) {
       if((rc = pthread_create(&threads[i], NULL, process_lines, (void *)(thread_data+i)))) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
       }
    }

    pthread_t writer_thread; //(pthread_t *)malloc(sizeof(pthread_t)); 
    //create writer thread
    if((rc = pthread_create(&writer_thread, NULL, write_results, (void *)(writer_data)))) {
         cout << "Error:unable to create thread," << rc << endl;
         exit(-1);
    }


    void *status;
    for(int i=0; i<num_threads; i++) {
      rc = pthread_join(threads[i], &status);
      if (rc) {
         printf("ERROR; return code from pthread_join() is %d\n", rc);
         exit(-1);
      }
     }

    rc = pthread_join(writer_thread, &status);
      if (rc) {
         printf("ERROR; return code from pthread_join() is %d\n", rc);
         exit(-1);
    }
}


void process_blastoutput(const Options& options, const GLOBAL_PARAMS &params) {
    unsigned int b =0;
    OutputParser parser(options, params);
    map<std::string, std::string> *annot_map = new map<std::string, std::string>;

    parser.initialize();

    THREAD_DATA *thread_data = new THREAD_DATA[options.num_threads];

    for(unsigned int i = 0; i < options.num_threads; i++) {
        thread_data[i].annot_map = annot_map;
        thread_data[i].options = options;
    }

    std::cout << "Read  annotation map \n";
    parser.create_annotation_dictionary(annot_map);
    
    std::cout << "Reading  refscores \n";
    parser.create_refBitScores(thread_data);

    
    std::cout << "Computing refscores \n";
    create_threads_refscores(options.num_threads, thread_data);


    WRITER_DATA *writer_data = new WRITER_DATA;

//    std::ofstream output(options.parsed_output.c_str(), std::ofstream::out);

    parser.initializeBatchReading();

    writer_data->output.open(options.parsed_output.c_str(), std::ofstream::out);

    writer_data->thread_data = thread_data;
    writer_data->num_threads = options.num_threads;

    while(parser.readABatch()) {  // main loop
       parser.distributeInput(thread_data);

       for(unsigned int i = 0; i < options.num_threads; i++) 
         thread_data[i].b = b;

       create_threads_parse(options.num_threads, thread_data, writer_data);
       b = (b+1)%2;
    } 

    for(unsigned int i = 0; i < options.num_threads; i++) 
       thread_data[i].b = b;
    write_results(writer_data);

    writer_data->output.close();
    parser.closeBatchReading();

}


void *write_results(void *_writer_data ) {

    std::cout << "Writing results \n";
    unsigned int b; 
    WRITER_DATA *writer_data = (WRITER_DATA *)_writer_data;
    unsigned int num_threads = writer_data->num_threads;
    THREAD_DATA *thread_data = writer_data->thread_data;
     

    vector< BLASTOUT_DATA >::iterator it;
    for(unsigned int i = 0; i < num_threads; i++) {
       b =  (thread_data[i].b+1)%2;
    //   std::cout << "Writing results from thread " << i << " buffer " << b << std::endl;
       for(it = thread_data[i].output[b].begin(); it != thread_data[i].output[b].end(); it++) {
           writer_data->output << it->query << "\t" <<  it->target << "\t" << it->q_length\
                  <<"\t" <<it->bitscore << "\t" << it->bsr <<"\t" << it->expect\
                  << "\t" << it->aln_length << "\t" << it->identity\
                  << "\t" << it->ec << "\t" << it->product <<  std::endl; 

       }
       thread_data[i].output[b].clear();
    }
    
    return (void *)NULL;

}
