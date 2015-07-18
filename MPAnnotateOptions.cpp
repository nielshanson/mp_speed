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
            this->blast_dir = argv[++i];
        }
        else if ( strncmp(argv[i], "--blast_dir", strlen("--blast_dir")) == 0 ) {
            this->blast_dir = argv[++i];
        }
        else if( strncmp(argv[i], "-a", strlen("-a")) == 0 ) {
            this->algorithm = argv[++i];
        }
        else if( strncmp(argv[i], "--algorithm", strlen("--algorithm")) == 0 ) {
            this->algorithm = argv[++i];
        }
        else if( strncmp(argv[i], "-f", strlen("-f")) == 0 ) {
            this->functional_categories = argv[++i];
        }
        else if( strncmp(argv[i], "--functional_categories", strlen("--functional_categories")) == 0 ) {
            this->functional_categories = argv[++i];
        }
        else if( strncmp(argv[i], "-i", strlen("-i")) == 0 ) {
            this->input_gff = argv[++i];
        }
        else if( strncmp(argv[i], "--input_gff", strlen("--input_gff")) == 0 ) {
            this->input_gff = argv[++i];
        }
        else if( strncmp(argv[i], "-s", strlen("-s")) == 0 ) {
            this->sample_name = argv[++i];
        }
        else if( strncmp(argv[i], "--sample_name", strlen("--sample_name")) == 0 ) {
            this->sample_name = argv[++i];
        }
        else if( strncmp(argv[i], "-r", strlen("-r")) == 0 ) {
            this->results_dir = argv[++i];
        }
        else if( strncmp(argv[i], "--results_dir", strlen("--results_dir")) == 0 ) {
            this->results_dir = argv[++i];
        }
        else if (strncmp(argv[i], "-p", strlen("-p")) == 0) {
            this->ptools_rxn_file = argv[++i];
        }
        else if (strncmp(argv[i], "--ptools_rxns", strlen("--ptools_rxns")) == 0) {
            this->ptools_rxn_file = argv[++i];
        }
        else if (strncmp(argv[i], "-t", strlen("-t")) == 0) {
            this->ptools_dir = argv[++i];
        }
        else if (strncmp(argv[i], "--ptools_dir", strlen("--ptools_dir")) == 0) {
            this->ptools_dir = argv[++i];
        }
        else if( strncmp(argv[i], "--tax", strlen("--tax")) == 0 ) {
            this->taxonomy = true;
        }
        else if( strncmp(argv[i], "-w", strlen("-w")) == 0 ) {
            this->weight_db = atof(argv[++i]);
        }
        else if( strncmp(argv[i], "--weight_db", strlen("--weight_db")) == 0 ) {
            this->weight_db = atof(argv[++i]);
        }
        else if( strncmp(argv[i], "--ncbi_catalog", strlen("--ncbi_catalog")) == 0 ) {
            this->ncbi_catalog = argv[++i];
        }
        else if( strncmp(argv[i], "--ncbi_catalog_names_map", strlen("--ncbi_catalog_names_map")) == 0 ) {
            this->ncbi_catalog_names_map = argv[++i];
        }
        else if ( strncmp(argv[i], "--ncbi_nodes", strlen("--ncbi_nodes")) == 0 ) {
            this->ncbi_nodes = argv[++i];
        }
        else if( strncmp(argv[i], "--num_threads", strlen("--num_threads")) == 0 ) {
            this->num_threads = atoi(argv[++i]);
        }
        else if( strncmp(argv[i], "--debug", strlen("--debug")) == 0 ) {
            this->debug = true;
        }
        else {
            cout << "ERROR: Cannot recognize argument " << argv[i] << endl;
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
         << "blast_dir:         " <<  blast_dir << "\n"\
         << "algorithm:         " <<  algorithm << "\n"\
         << "functional_categories: " <<  functional_categories << "\n"\
         << "input_gff:         " << input_gff << "\n"\
         << "sample_name:       " <<  sample_name << "\n"\
         << "taxonomy:          " <<  taxonomy << "\n"\
         << "results_dir:       " <<  results_dir << "\n"\
         << "ptools_rxn_file:   " <<  ptools_rxn_file << "\n"\
         << "ptools_dir:        " <<  ptools_dir << "\n"\
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
