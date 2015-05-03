//
// Created by Niels Hanson on 2015-04-29.
//

/*
 * annotation.c++: This class library contains functions associated the processing of multiple parsed database results.
 * (MPAnnotate)
 */

#include "annotation.h"

using namespace std;

#define PRINT_INTERVAL 10000

void readContigLengths(string file, map<string, unsigned int> &contig_lengths) {


    string filename = file;
    std::ifstream input;
    char buf[1000];
    vector <char *> fields;
    int count = 0;

    cout << "Reading refscores " <<  endl;
    cout << "Filename " << filename << "\n";
    input.open(filename.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
        return ;
    }

    // int *counts  = (int *)calloc(options.num_threads, sizeof(int));
    string line;
    string contig_id;

    while( std::getline(input, line ).good()) {
        split(line, fields, buf,'\t');
        if( fields.size()!=3) {
            contig_id.clear();
            input.close();
            return;
        }

        contig_id = fields[0];
        contig_lengths[contig_id] = atoi(fields[2]);



        if (count%PRINT_INTERVAL==0)
            std::cout << "x " << count << std::endl;
        count++;
    }
    input.close();



    std::cout << "Number of contig lengths loaded " <<  count << std::endl;


}


/*
 * Given a product string this function computes the products word information score. I.e., one point for every
 * non-trivial descripive word.
 */
float wordInformation(string product) {

    // Split product into individual words
    char buf[10000]; // TODO: Check to see if this buffer is needed.
    vector <char *> words;
    split(product, words, buf,' ');

    float word_information_score = 0.0; // information score of the current words

    // Prepare map of stop words to check for
    static const string arr[] = {"", "is", "have", "has", "will", "can", "should",  "in",
                                 "at", "upon", "the", "a", "an", "on", "for", "of", "by",
                                 "with" ,"and",  ">", "predicted", "protein", "conserved" };
    map <string, int> stop_words;
    for (int i = 0; i < arr->size(); i++) {
        stop_words[arr[i]] = 1;
    }

    // Check each word for membership in stop_words and presence of underscores '_'
    for (int i = 0; i < words.size(); i++) {
        string my_word = words[i];
        if (stop_words.count(my_word) <= 0) {
            // not a stop word
            if (my_word.find("_") == std::string::npos) {
                // does not contain an underscore
                word_information_score += 1.0;
            }
        }
    }

    return word_information_score;
}


/*
 * Computes the annotation value of the given ANNOTATION object. Returns an annotation value based on the presence
 * of Enzyme Commission (EC) numbers (+10) and the number of non-trivial words in the ANNOTATION product field.
 */
float computeAnnotationValue(ANNOTATION annotation) {

    float score = 0.0; // overall annotation score

    if (annotation.ec != "") {
        // annotation object has non-trivial EC number
        score += 10;
    }

    if (annotation.product.find("hypothetical protein") != std::string::npos) {
        score += wordInformation(annotation.product);
    }

    return score;

}

int processParsedBlastout(string db_name, float weight, string blastoutput, MPAnnotateOptions options, map<string, ANNOTATION> &annotation_results) {

    if (options.debug) cout << "In processParsedBlastout()" << endl;

    // Prepare inputs and buffers
    string filename = options.blast_dir + "/" + blastoutput; // BLAST/LASTout.parsed.txt
    std::ifstream input; // Input filestream
    char buf[10000]; // Temp buffer TODO check to see if multiple buffers nessisary
    vector <char *> fields; // Vector for parsed fields
    int count = 0; // Line count
    char tempbuf[1000];

    if (options.debug) {
        cout << "Reading ParsedBlastout: " << filename << "\n";
    }

    // Open the file.
    input.open(filename.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '" << filename << "'. Bailing out." << std::endl;
        exit(1);
    }

    // ANNOTATION fields for parsing BLAST/LASTout.annotated.txt files
    ANNOTATION annotation;
    string line;
    string query_id;
    string bsr;
    string ec;
    string product;
    string taxonomy = "";
    float value;

    // Parse each line and create ANNOTATION objects
    while( std::getline(input, line ).good()) {
        split(line, fields, buf,'\t');
        if (count == 0) { count++; continue; };

        if( fields.size()!=10) {
            // Not a parsed blast file
            cerr << "Parsed BLAST/LASTout file " <<  filename << " did not have the 10 columns." << endl;
            query_id.clear();
            input.close();
            return 1;
        }

        query_id = fields[0];
        bsr = fields[4];
        ec = fields[8];
        product = fields[9];

        int i =0;
        for(i =0; i < db_name.size(); i++ )
            tempbuf[i] = std::toupper(db_name[i]);
        tempbuf[i]='\0';
        db_name = string(tempbuf);


        // Construct annotation
        ANNOTATION annotation;
        annotation.bsr = std::stof(bsr);
        annotation.ec = ec;
        annotation.product = product;

        // Extract RefSeq taxonomy from product field
        if (options.taxonomy && (db_name.find("REFSEQ") != std::string::npos)) {
            // TODO: Could be more reliable to get taxonomy from GI number.
            taxonomy = getTaxonomyFromProduct(product.c_str());
        }
        annotation.taxonomy = taxonomy;


        // Compute information score of current annotation.
        annotation.value = computeAnnotationValue(annotation) * weight;

        // Check to see if current query_id is already in this database's annotation_results
        if (annotation_results.count(query_id) <= 0) {
            annotation_results[query_id] = annotation;
        }
        else {
            // Replace existing ANNOTATION with the current annotation if strong annotation score
            if (annotation_results[query_id].value < annotation.value) {
                annotation_results[query_id] = annotation;
            }
        }

        if (options.debug) {
            if (count%PRINT_INTERVAL==0)
                std::cout << "x " << count << std::endl;
        }

        count++;
    }
    input.close();

    if (options.debug) {
        std::cout << "Number of parsed BLAST/LAST results loaded " <<  count << std::endl;
    }
    return count;
}

/*
 * Given a string with the location of the blast_results directory, getBLASTFileNames returns
 * a list of the out.parsed.txt files.
 */
int getBlastFileNames(string blastdir, string sample_name, MPAnnotateOptions options, DB_INFO &db_info) {

    regex_t regex; // Regular expression.
    char tempbuf[1000]; // Buffer.

    // Convert algorithm name to uppercase (i.e., BLAST, LAST)
    string algorithm = options.algorithm;
    int i =0;
    for(i =0; i < algorithm.size(); i++ )
        tempbuf[i] = std::toupper(algorithm[i]);
    tempbuf[i]='\0';
    algorithm = string(tempbuf);

    // Create regular expression pattern for BLAST/LASTout.parsed.txt files
    string regexp_text = sample_name + "[.](.*)[.]" + algorithm + "out.parsed.txt";
    compile_regex(&regex, regexp_text.c_str());

    // Open directory
    vector<string> files; // list of files in directory
    DIR *dp; // directory pointer
    struct dirent *dirp;
    if((dp  = opendir(blastdir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening cirectory " << blastdir << endl;
        return errno;
    }

    // Add files in directory to files vector
    while ((dirp = readdir(dp)) != NULL) {
        files.push_back(string(dirp->d_name));
    }
    closedir(dp);

    // Check to see if files match parsed.txt pattern
    string database;
    for (unsigned int i = 0;i < files.size();i++) {
        database = getpattern(&regex, files[i].c_str(), 1) ;
        if( database.size() > 0)  {
            // Pattern found: add database name, filename, and weight to db_info
            db_info.db_names.push_back(database);
            db_info.input_blastouts.push_back(files[i]);
            db_info.weight_dbs.push_back(1.0);
        }
    }

    return 0;
}

/*
 *
 */
void createAnnotation(map<string, float> dbname_weight, ANNOTATION_RESULTS results_dictionary, MPAnnotateOptions options, map<string, unsigned int> contig_lengths) {
    // create_annotation(dbname_weight, results_dictionary, opts.input_gff, opts.rRNA_16S, opts.tRNA, opts.output_gff, opts.output_comparative_annotation, contig_lengths)
    // orf_dictionary={};

    if (options.debug) {
        cout << "In createAnnotation()" << endl;
    }
    string input_gff = options.input_gff;
    string rRNA_16S = options.rRNA_16S;
    string tRNA = options.tRNA;
    string output_gff = options.output_gff;
    string output_comp_annot = options.output_comp_annot;

    // read input ORFs the input GFF file
    cout << options.input_gff << endl;

    // Prepare inputs and buffers
    string filename = options.blast_dir + "/" + blastoutput; // BLAST/LASTout.parsed.txt
    std::ifstream input; // Input filestream
    char buf[10000]; // Temp buffer TODO check to see if multiple buffers nessisary
    vector <char *> fields; // Vector for parsed fields
    int count = 0; // Line count
    char tempbuf[1000];


//    output_gff_tmp = output_gff + ".tmp";
//    outputgff_file = open( output_gff_tmp, 'w');
//    output_comp_annot_file1 = open( output_comparative_annotation + '.1.txt', 'w');
//    output_comp_annot_file2 = open( output_comparative_annotation + '.2.txt', 'w');
//
//    output_comp_annot_file1_Str = 'orf_id\tref dbname\tEC\tproduct\tvalue';
//    fprintf(output_comp_annot_file1,'%s\n', output_comp_annot_file1_Str);
//
//    output_comp_annot_file2_Str = 'orf_id'
//    dbnames = dbname_weight.keys()
//    for dbname in dbnames:
//    weight = dbname_weight[dbname]
//    output_comp_annot_file2_Str += '\t{0}(EC) \t{0}(product)\t{0}(value)'.format(dbname)
//    fprintf(output_comp_annot_file2,'%s\n', output_comp_annot_file2_Str)
}