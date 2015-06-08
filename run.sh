#!/usr/bin/env bash
#parse -d metacyc-v4-2011-07-03 -b /Volumes/3TB/SPEED/output/mgm4541646/blast_results//mgm4541646.metacyc-v4-2011-07-03.LASTout -m /Users/sgeadmin/MetaPathways/databases//functional/formatted/metacyc-v4-2011-07-03-names.txt -r /Volumes/3TB/SPEED/output/mgm4541646/blast_results//mgm4541646.refscores.LAST --min_bsr 0.4 --min_score 20 --min_length 30 --max_evalue 0.000001 -a LAST --num_threads 20 -o x.txt
#parse -d metacyc-v4-2011-07-03 -b /Volumes/3TB/SPEED/output/mgm4541646/blast_results//mgm4541646.metacyc-v4-2011-07-03.LASTout -m /Users/sgeadmin/MetaPathways/databases//functional/formatted/seed-2014-01-30-names.txt -r /Volumes/3TB/SPEED/output/mgm4541646/blast_results//mgm4541646.refscores.LAST --min_bsr 0.4 --min_score 20 --min_length 60 --max_evalue 0.000001 -a LAST --num_threads 20
#parse -d seed-2014-01-30 -b /Volumes/3TB/PROCESSED_DATA/DATA/SHORT_READS/mgm4510174.seed-2014-01-30.LASTout -m /Users/sgeadmin/MetaPathways/databases//functional/formatted/seed-2014-01-30-names.txt -r /Volumes/4TB/PROCESSED_DATA/DATA/SHORT_READS/output/mgm4510174/blast_results/mgm4510174.refscores.LAST  --min_bsr 0.4 --min_score 20 --min_length 60 --max_evalue 0.000001 -a LAST --num_threads 20 -o x.txt

#parse -d metacyc-v4-2011-07-03 -b  /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/blast_results/hmp_airways_SRS014682.metacyc-v4-2011-07-03.LASTout -m /Users/sgeadmin/MetaPathways/databases/functional/formatted/metacyc-v4-2011-07-03-names.txt -r /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/blast_results/hmp_airways_SRS014682.refscores.LAST --min_bsr 0.4 --min_score 20 --min_length 30 --max_evalue 0.000001 --algorithm LAST --num_threads 1 -o x.txt -s hmp_airways_SRS014682


# kishori's original
#mp_annotate --input_gff /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/orf_prediction/hmp_airways_SRS014682.unannot.gff \
# --output_gff /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/results/annotation_table/hmp_airways_SRS014682 \
# --output_comparative_annotation /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/results/annotation_table/hmp_airways_SRS014682 \
# -D /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/blast_results/ \
# -s hmp_airways_SRS014682 \
# -m /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/preprocessed/hmp_airways_SRS014682.mapping.txt \
# -a LAST \
# --tax \
# --debug \
# --num_threads 1

# ./mp_annotate --input_gff data/mp_output/hmp_airways_SRS014682/orf_prediction/hmp_airways_SRS014682.unannot.gff \
# --output_gff data/mp_output/hmp_airways_SRS014682/results/annotation_table/hmp_airways_SRS014682 \
# --output_comparative_annotation data/mp_output/hmp_airways_SRS014682/results/annotation_table/hmp_airways_SRS014682 \
# -D data/mp_output/hmp_airways_SRS014682/blast_results/ \
# -s hmp_airways_SRS014682 \
# -m data/mp_output/hmp_airways_SRS014682/preprocessed/hmp_airways_SRS014682.mapping.txt \
# -a LAST \
# --tax \
# --debug \
# --num_threads 1


# Test for create_ptools_input
time ./mp_create_ptools_input --ptools_rxns /Users/nielshanson/Dropbox/projects/mp_speed/data/metacyc_enzymes_rxns_ecs.txt \
                         --anno_table /Users/nielshanson/Dropbox/projects/mp_speed/data/hmp_airways_SRS014682.functional_and_taxonomic_table.txt3 \
                         --ptools_dir /Users/nielshanson/Dropbox/projects/mp_speed/data


#for f in  ALCAME_VIR_PRO  CARD_ABR CARD_Anti_Biotic_Resis CAZY_2014_09_04 COG_2013-12-27  kegg-pep-2011-06-18  refseq-nr-2014-01-18  seed-2014-01-30; 
# do 
#   echo $f
#parse -d ${f} -b  /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/blast_results/hmp_airways_SRS014682.${f}.LASTout -m /Users/sgeadmin/MetaPathways/databases/functional/formatted/${f}-names.txt -r /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/blast_results/hmp_airways_SRS014682.refscores.LAST --min_bsr 0.4 --min_score 20 --min_length 30 --max_evalue 0.000001 --algorithm LAST --num_threads 10 -o /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/blast_results/hmp_airways_SRS014682.${f}.LASTout.parsed.txt -s hmp_airways_SRS014682
#done;

#parse -d metacyc-v4-2011-07-03 -b  /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/blast_results/hmp_airways_SRS014682.metacyc-v4-2011-07-03.LASTout -m /Users/sgeadmin/MetaPathways/databases/functional/formatted/metacyc-v4-2011-07-03-names.txt -r /Volumes/3TB/SPEED/output/hmp_airways_SRS014682/blast_results/hmp_airways_SRS014682.refscores.LAST --min_bsr 0.4 --min_score 20 --min_length 30 --max_evalue 0.000001 --algorithm LAST --num_threads 1 -o x.txt -s hmp_airways_SRS014682

