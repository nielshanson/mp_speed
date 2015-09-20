#!/usr/bin/env bash

# integration_run.sh
# Script to test all C++ code in MetaPathways v3.0

#data_folder=/home/ubuntu/workspace/data
data_folder=/Users/nielshanson/Dropbox/projects/mp_speed/data
mp_output=${data_folder}/mp_output
mp_databases=${data_folder}/mp_databases/MetaPathways_DBs
mp_source=${data_folder}/metapathways2
resources=${data_folder}/resources

# sample paramters
sample_name=hmp_airways_SRS014682
num_threads=4

## Remove old data

# # LASTDBs
# my_cmd="rm ${mp_databases}/functional/formatted/*.bck"
# echo $my_cmd
# eval $my_cmd
# my_cmd="rm ${mp_databases}/functional/formatted/*.des"
# echo $my_cmd
# eval $my_cmd
# my_cmd="rm ${mp_databases}/functional/formatted/*.prj"
# echo $my_cmd
# eval $my_cmd
# my_cmd="rm ${mp_databases}/functional/formatted/*.sds"
# echo $my_cmd
# eval $my_cmd
# my_cmd="rm ${mp_databases}/functional/formatted/*.ssp"
# echo $my_cmd
# eval $my_cmd
# my_cmd="rm ${mp_databases}/functional/formatted/*.suf"
# echo $my_cmd
# eval $my_cmd
# my_cmd="rm ${mp_databases}/functional/formatted/*.tis"
# echo $my_cmd
# eval $my_cmd

# # LASTouts
# my_cmd="rm ${mp_output}/${sample_name}/blast_results/*LASTout*"
# echo $my_cmd
# eval $my_cmd

# # functional hierarchy hits
# my_cmd="rm ${mp_output}/${sample_name}/results/annotation_table/*tree.count.txt"
# echo $my_cmd
# eval $my_cmd

# # ptools/
# my_cmd="rm ${mp_output}/${sample_name}/ptools/*"
# echo $my_cmd
# eval $my_cmd


# ## Integration run
# # Format lastdbs
# for db_file in ${mp_databases}/functional/*-*
# do
#     filename="${db_file##*/}"
#     my_cmd="./lastdb+ -p ${mp_databases}/functional/formatted/${filename} $db_file"
#     echo $my_cmd
#     eval $my_cmd
# done

# # LAST/data
# for db_file in ${mp_databases}/functional/*-*
# do
#     filename="${db_file##*/}"
#     my_cmd="./lastal+ -S 20 \
#                       -E 1e-05 \
#                       -P ${num_threads} \
#                       -K 10 \
#                       -f 2 \
#                       -o ${mp_output}/${sample_name}/blast_results/${sample_name}.${filename}.LASTout \
#                       ${mp_databases}/functional/formatted/${filename} \
#                       ${mp_output}/${sample_name}/orf_prediction/${sample_name}.qced.faa"
#     echo $my_cmd
#     eval $my_cmd
# done


# parseB/LAST
# for file in ${mp_output}/${sample_name}/blast_results/*LASTout
# do
#     # extract database names from LASTout file
#     db="$(echo ${file} | sed -e 's/\(.*\)\.LASTout$/\1/g' | sed -e 's/.*\.//g')"
#
#     # Parse B/LAST output
#     my_cmd="./mp_parseblast -d ${db} \
#           -b ${mp_output}/${sample_name}/blast_results/${sample_name}.${db}.LASTout \
#           -m ${mp_databases}/functional/formatted/${db}-names.txt \
#           -r ${mp_output}/${sample_name}/blast_results/${sample_name}.refscores.LAST \
#           --min_bsr 0.4 \
#           --min_score 20 \
#           --min_length 30 \
#           --max_evalue 0.000001 \
#           -a LAST \
#           --num_threads ${num_threads} \
#           -o ${mp_output}/${sample_name}/blast_results/${sample_name}.${db}.LASTout.parsed.txt"
#     echo $my_cmd
#     eval $my_cmd
# done

# Annotate
#my_cmd="./mp_annotate --input_gff ${mp_output}/${sample_name}/orf_prediction/${sample_name}.unannot.gff \
#                    --results_dir ${mp_output}/${sample_name}/results/annotation_table/ \
#                    --blast_dir ${mp_output}/${sample_name}/blast_results/ \
#                    --functional_categories ${mp_databases}/functional_categories \
#                    --ptools_rxns ${resources}/metacyc_enzymes_rxns_ecs.txt \
#                    --ptools_dir ${mp_output}/${sample_name}/ptools \
#                    --sample_name ${sample_name} \
#                    --algorithm LAST \
#                    --tax \
#                    --ncbi_catalog_map ${resources}/RefSeq-release69.catalog.small.txt \
#                    --ncbi_catalog_names_map ${resources}/RefSeq-release69.catalog.taxid2taxa.txt \
#                    --ncbi_nodes ${resources}/ncbi_nodes_parent_child_ids.txt \
#                    --debug \
#                    --num_threads ${num_threads}"
#echo $my_cmd
#eval $my_cmd

my_cmd="./mp_annotate --input_gff /Users/nielshanson/Documents/pycharm/metapathways3/regtests/mp3_benchmarking/mp3_output/bench_90k/orf_prediction//bench_90k.unannot.gff \
 --results_dir /Users/nielshanson/Documents/pycharm/metapathways3/regtests/mp3_benchmarking/mp3_output/bench_90k/results//annotation_table// \
 --blast_dir /Users/nielshanson/Documents/pycharm/metapathways3/regtests/mp3_benchmarking/mp3_output/bench_90k/blast_results/ \
 --functional_categories /Users/nielshanson/Documents/pycharm/metapathways3/databases//functional_categories \
 --ptools_rxns /Users/nielshanson/Documents/pycharm/metapathways3/databases//functional_categories/metacyc_enzymes_rxns_ecs.txt \
 --ptools_dir /Users/nielshanson/Documents/pycharm/metapathways3/regtests/mp3_benchmarking/mp3_output/bench_90k/ptools/ \
 --sample_name bench_90k \
 --algorithm LAST \
 --tax \
 --ncbi_catalog_map /Users/nielshanson/Documents/pycharm/metapathways3//resources//RefSeq-release69.catalog.small.txt \
 --ncbi_catalog_names_map /Users/nielshanson/Documents/pycharm/metapathways3//resources//RefSeq-release69.catalog.taxid2taxa.txt \
 --ncbi_nodes /Users/nielshanson/Documents/pycharm/metapathways3//resources//ncbi_nodes_parent_child_ids.txt \
 --debug \
 --num_threads 4"

echo $my_cmd
eval $my_cmd

#my_cmd="./mp_annotate --input_gff /Users/nielshanson/Desktop/kidney_stones_out/ALaw-C5_S6_contigs/orf_prediction/ALaw-C5_S6_contigs.unannot.gff \
# --results_dir /Users/nielshanson/Desktop/kidney_stones_out/ALaw-C5_S6_contigs/results/annotation_table/ \
# --blast_dir /Users/nielshanson/Desktop/kidney_stones_out/ALaw-C5_S6_contigs/blast_results/ \
# --functional_categories /Users/nielshanson/Documents/pycharm/metapathways3/databases/functional_categories \
# --ptools_rxns /Users/nielshanson/Documents/pycharm/metapathways3/databases/functional_categories/metacyc_enzymes_rxns_ecs.txt \
# --ptools_dir /Users/nielshanson/Desktop/kidney_stones_out/ALaw-C5_S6_contigs/ptools/ \
# --sample_name ALaw-C5_S6_contigs \
# --algorithm LAST \
# --tax \
# --ncbi_catalog_map /Users/nielshanson/Documents/pycharm/metapathways3/resources/RefSeq-release69.catalog.small.txt \
# --ncbi_catalog_names_map /Users/nielshanson/Documents/pycharm/metapathways3/resources/RefSeq-release69.catalog.taxid2taxa.txt \
# --ncbi_nodes /Users/nielshanson/Documents/pycharm/metapathways3/resources/ncbi_nodes_parent_child_ids.txt \
# --debug \
# --num_threads 4"
#
#echo $my_cmd
#eval $my_cmd


# Prepare Ptools Input
#./mp_create_ptools_input --ptools_rxns /Users/nielshanson/Dropbox/projects/mp_speed/data/metacyc_enzymes_rxns_ecs.txt \
#                         --anno_table /Users/nielshanson/Dropbox/projects/mp_speed/data/hmp_airways_SRS014682.functional_and_taxonomic_table_copy.txt \
#                         --ptools_dir /Users/nielshanson/Dropbox/projects/mp_speed/data