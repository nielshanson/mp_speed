//
// Created by Niels Hanson on 2015-04-29.
//

#include "MPAnnotateOptions.h"


/*
 * Parse arguments for MetaPathways annotation step.
 */
bool MPAnnotateOptions::SetOptions(int argc, char *argv[]) {
    for(int i = 1; i < argc ; i++) {
        if( strncmp(argv[i], "-b", strlen("-b")) == 0 ) {
            this->input_blastout = argv[++i];
        }
        else if( strncmp(argv[i], "-a", strlen("-a")) == 0 ) {
            this->algorithm = argv[++i];
        }
        else if( strncmp(argv[i], "-m", strlen("-m")) == 0 ) {
            this->contig_map_file = argv[++i];
        }
        else if(strncmp(argv[i], "-r", strlen("-r")) == 0 ) {
            this->database_name =argv[++i];
        }
        else if(strncmp(argv[i], "-D", strlen("-D")) == 0 ) {
            this->blast_dir =argv[++i];
        }
        else if( strncmp(argv[i], "--input_gff", strlen("--input_gff")) == 0 ) {
            this->input_gff = argv[++i];
        }
        else if( strncmp(argv[i], "--rRNA_16S", strlen("--rRNA_16S")) == 0 ) {
            this->rRNA_16S = argv[++i];
        }
        else if( strncmp(argv[i], "--tRNA", strlen("--tRNA")) == 0 ) {
            this->tRNA = argv[++i];
        }
        else if( strncmp(argv[i], "--tax", strlen("--tax")) == 0 ) {
            this->taxonomy = true;
        }
        else if( strncmp(argv[i], "--output_gff", strlen("--output_gff")) == 0 ) {
            this->output_gff = argv[++i];
        }
        else if( strncmp(argv[i], "--output_comparative_annotation", strlen("--output_comparative_annotation")) == 0 ) {
            this->output_comp_annot = argv[++i];
        }
        else if( strncmp(argv[i], "-s", strlen("-s")) == 0 ) {
            this->sample_name = argv[++i];
        }
        else if( strncmp(argv[i], "-s", strlen("-s")) == 0 ) {
            this->weight_db = atof(argv[++i]);
        }
        else if( strncmp(argv[i], "--num_threads", strlen("--num_threads")) == 0 ) {
            this->input_gff = atoi(argv[++i]);
        }
        else if( strncmp(argv[i], "--debug", strlen("--debug")) == 0 ) {
            this->debug = true;
        }
        else {
            cout << "ERROR: Cannot recognize argument " << argv[i] << std::endl;;
            return false;
        }
    } //for loop for arguments processing

    return true;

};


/*
 * Print current options
 */

void MPAnnotateOptions::printOptions() {
    cout << "CURRENT OPTIONS:\n"\
         << "---------------\n"
         << "input_blastout:    " <<  input_blastout << "\n"\
         << "algorithm:         " <<  algorithm << "\n"\
         << "contig_map_file:   " <<  contig_map_file << "\n"\
         << "database_name:     " <<  database_name << "\n"\
         << "blast_dir:         " <<  blast_dir << "\n"\
         << "rRNA_16S:          " <<  rRNA_16S  << "\n"\
         << "tRNA:              " <<  tRNA << "\n"\
         << "sample_name:       " <<  sample_name << "\n"\
         << "taxonomy:          " <<  taxonomy << "\n"\
         << "output_gff:        " <<  output_gff << "\n"\
         << "output_comp_annot: " <<  output_comp_annot << "\n"\
         << "input_gff:         " <<  input_gff << "\n"\
         << "num_threads:       " <<  num_threads << "\n"\
         << "weight_db:         " <<  weight_db << "\n"\
         << "debug:             " <<  debug << "\n"
         << "---------------\n\n";
    
}

/*
 * Print usage information.
 */
void MPAnnotateOptions::printUsage(char *arg) {
    std::cout << "USAGE : "   << arg <<"\n"\
              << "      : -b  <blastoutput>                                               [REQUIRED]\n"\
              << "      : -d  <database name>                                             [REQUIRED]\n"\
              << "      : -m  <contig_map_file>                                           [REQUIRED]\n"\
              << "      : -r  <refscorefile>                                              [REQUIRED]\n"
              << "      : -D  <blastdir>                                                  [REQUIRED]\n"\
              << "      : --output_comparative_annotation <comparative_annotation_table>  [REQUIRED]\n"\
              << "      : --input_gff <input_gff_file>                                    [REQUIRED] \n"\
              << "      : -s <sample_name>                                                [REQUIRED] \n"\
              << "      : -a  <algorithm>                                 [OPTIONAL, default: BLAST]\n"\
              << "      : --rRNA_16S <16s_rRNA_stats_file>                [OPTIONAL]\n"\
              << "      : --tRNA <tRNA_rRNA_stats_file>                   [OPTIONAL]\n"\
              << "      : --tax                                           [OPTIONAL, default true]\n"\
              << "      : --output_gff <output_gff_file>                  [OPTIONAL, default 30]\n"\
              << "      : --num_threads <t>                               [OPTIONAL, default 1]\n"\
              << std::endl;
};
