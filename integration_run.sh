#!/usr/bin/env bash

# integration_run.sh
# Script to test all C++ code in MetaPathways v3.0

data_folder=/Users/nielshanson/Dropbox/projects/mp_speed/data
mp_output=${data_folder}/mp_output
mp_databases=${data_folder}/mp_databases/MetaPathways_DBs
mp_source=${data_folder}/metapathways2
resources=${data_folder}/resources

# sample paramters
sample_name=hmp_airways_SRS014682
num_threads=4

# LAST

# Parse B/LAST output
my_cmd="./parse -d metacyc-v4-2011-07-03 \
      -b ${mp_output}/${sample_name}/blast_results/${sample_name}.metacyc-v4-2011-07-03.LASTout \
      -m ${mp_databases}/functional/formatted/metacyc-v4-2011-07-03-names.txt \
      -r ${mp_output}/${sample_name}/blast_results/${sample_name}.refscores.LAST \
      --min_bsr 0.4 \
      --min_score 20 \
      --min_length 30 \
      --max_evalue 0.000001 \
      -a LAST \
      --num_threads ${num_threads} \
      -o ${mp_output}/${sample_name}/blast_results/${sample_name}.metacyc-v4-2011-07-03.LASTout.parsed.txt"

echo $my_cmd
# eval $my_cmd

# Annotate
my_cmd="./mp_annotate --input_gff ${mp_output}/${sample_name}/orf_prediction/${sample_name}.unannot.gff \
                    --output_gff ${mp_output}/${sample_name}/orf_prediction/${sample_name}.gff \
                    --output_comparative_annotation ${mp_output}/${sample_name}/results/annotation_table/${sample_name}.functional_and_taxonomic_table.txt \
                    -D ${mp_output}/${sample_name}/blast_results/ \
                    -s ${sample_name} \
                    -m ${mp_output}/${sample_name}/preprocessed/${sample_name}.mapping.txt \
                    -a LAST \
                    --tax \
                    --debug \
                    --num_threads ${num_threads}"

#echo $my_cmd
#eval $my_cmd

# Prepare Ptools Input
#./mp_create_ptools_input --ptools_rxns /Users/nielshanson/Dropbox/projects/mp_speed/data/metacyc_enzymes_rxns_ecs.txt \
#                         --anno_table /Users/nielshanson/Dropbox/projects/mp_speed/data/hmp_airways_SRS014682.functional_and_taxonomic_table_copy.txt \
#                         --ptools_dir /Users/nielshanson/Dropbox/projects/mp_speed/data