//
// Created by Niels Hanson on 2015-05-12.
//

#include "MPCreatePToolsInput.h"

using namespace std;

int main( int argc, char** argv) {

    // Parse options
    MPCreatePToolsInputOptions options;
    if ( ! options.SetOptions(argc, argv) ) {
        // options.printUsage(argv[0]);
        exit(1);
    }

    vector<string> ptools_list;
    processPToolsRxnsFile(options.ptools_rxn_file, ptools_list);

    PTOOLS_NODE *root = new PTOOLS_NODE();
    PTOOLS_NODE *start;

    // Iterate through the ptools list of annotations
    char buf[10000]; // Temp buffer
    std::vector<char *> words; // Vector for parsed fields

    // Parse each line of ptools annotations
    for (std::vector<string>::iterator itr = ptools_list.begin() ; itr != ptools_list.end(); ++itr) {
        split(*itr, words, buf,' ');
        start = root;
        for( unsigned int i=0; i < words.size(); i++) {
            string word = string(words[i]);
            // cout << word << endl;
            if(!start->hasChild(word)) {
                // Create child node
                start->insertChild(word);
                if (i == (words.size()-1)) {
                    // Flag node as finished
                    start->flagChildFinished(word);
                }
            }
            start = start->getChildNode(word);
        }
    }
    cout << "Loaded tree" << endl;
    //printMetaCycTree(root);

    LIST *my_list = new LIST(root); // Linked list to store MetaCyc tree pointers
    processAnnotationsForPTools(my_list, root, options.annotation_table, options.ptools_dir);

    exit(1);

//    insert into list
//    my_list->insert("One",1);
//    my_list->insert("Two",1);
//    my_list->insert("Three",1);
//
//    // Move current to next node
//    my_list->nextNode();
//    // Test delete
//    my_list->deleteCurr();
//    // Test insert
//    my_list->insertAtCurr("Four",1);
//
//    // Print out list
//    LIST_NODE *list_itr = my_list->head;
//    while(list_itr != NULL) {
//        cout << list_itr->my_data << endl;
//        list_itr = list_itr->next;
//    }


    // print out ptools tree
    // start = head;
    // print_dfs(start, "");

    // Create listnode linked list


}

/*
 * Print out MetaCyc word trie using a stack-based recursive post-order traversal
 */
void printMetaCycTree(PTOOLS_NODE *root) {
    vector<PTOOLS_NODE *> node_stack; // Stack of nodes
    vector<string> path; // Stack for path tracing
    string headpath = ""; // Path string

    // Push first nodes
    node_stack.push_back(root);
    path.push_back(headpath);

    while(!node_stack.empty()) {
        // Pop stacks
        PTOOLS_NODE* top = node_stack.back();
        node_stack.pop_back();
        headpath = path.back(); // Get current path
        path.pop_back();

        if (top->complete==true) {
            // Complete MetaCyc annotation
            cout << headpath << endl;
        }
        // Push children to stacks
        map<string, PTOOLS_NODE*> children = top->children;
        for (std::map<string, PTOOLS_NODE*>::iterator itr=children.begin(); itr != children.end(); ++itr) {
            node_stack.push_back(itr->second);
            path.push_back(headpath + itr->first + " "); // Keep track of path
        }
    }
}

/*
 * Prints out trie in using depth-first-search post order traversal.
 */
void print_dfs(PTOOLS_NODE *node, string line) {
    if (node->children.size() == 0) {
        // at leaf
        cout << line << endl;
    } else {
        for (std::map<string, PTOOLS_NODE*>::iterator itr=node->children.begin(); itr != node->children.end(); ++itr) {
            print_dfs(itr->second, line + " " + itr->first);
        }
    }
}

/*
 * Scans Pathway Tools reaction file and puts reaction strings into ptools_list
 */
void processPToolsRxnsFile( string ptools_rxn_file, vector<string> &ptools_list ){
    string filename  = ptools_rxn_file;
    std::cout << "Reading ptools_rxn_file " << filename <<  std::endl;
    input.open(filename.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '"<< filename <<"'. Bailing out." << std::endl;
        return ;
    }

    // Scanning variables
    int count =0;
    string line;
    char buf[10000]; // Temp buffer
    vector <char *> fields; // Vector for parsed fields

    while( std::getline(input, line ).good()) {
        // Skip header
        if (count == 0) {
            count++;
            continue;
        }

        split(line, fields, buf,'\t');
        ptools_list.push_back(fields[1]); // Ptools enzyme description field.
        count++;
    }

    input.close();
}

/*
 * Scans annotations and checks to see if annotations (or their sub-strings) are complete annotations
 * in the MetaCyc trie.
 */
void processAnnotationsForPTools(LIST *my_list, PTOOLS_NODE *root, string annotation_file, string ptools_dir) {

    // Opening annotation and output file
    std::cout << "Reading annotation_file " << annotation_file <<  std::endl;
    input.open(annotation_file.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '"<< annotation_file <<"'. Bailing out." << std::endl;
        return ;
    }
    string pf_file = ptools_dir + "/" + "0.pf";
    std::cout << "Opening handle to .pf file " << pf_file <<  std::endl;

    output.open(pf_file.c_str(), std::ofstream::out);
    if(!output.good()){
        std::cerr << "Error opening '"<< pf_file <<"'. Bailing out." << std::endl;
        return ;
    }

    // Scanning variables
    int count =0;
    string line;
    char buf[10000]; // Temp buffer
    vector <char *> fields; // Vector for parsed fields
    vector <char *> annotation_words; // Vector of annotation words

    // .pf file fields
    string orf_id; // MetaPathways ORFID
    string annotation_product; // MetaCyc FUNCTION field
    string ec_number; // MetaCyc EC field
    int start_base = 0; // Starting base
    int length = 0; // Basenumber

    // string metacyc_rxn_id; // MetaCyc reaction frame ID
    vector<string> annotations;

    // annotation variables
    PTOOLS_NODE *ptools_ptr = root;
    bool complete = false;
    vector <string> word_list;
    vector <string> max_word_list;
    string annotation = "";

    // For each annotation in annotation_file
    while( std::getline( input, line ).good() ) {
        if (count == 0) {
            count++;
            continue; // skip header
        }

        // Clear annotation fields
        orf_id.clear();
        annotation_product.clear();
        ec_number.clear();
        length = -1;
        fields.clear();

        // String line into fields and split annotation into words
        split(line, fields, buf, '\t');

        // Fill annotations
        orf_id = fields[0];
        annotation_product = fields[9];
        ec_number = fields[7];
        length = atoi(fields[1]);

        // Split annotation into separate words
        split(fields[9], annotation_words, buf, ' ');

        // Prepare annotation constants
        ptools_ptr = root;
        complete = false;
        word_list.clear();
        max_word_list.clear();
        annotation.clear();

        // Process annotation through MetaCyc trie
        // annotation_product = processAnnotationForPtools(annotation_words, root);
        annotation_product = processAnnotationForPtools(annotation_words, root, ptools_ptr, complete, word_list, max_word_list, annotation);

        // If valid annotation product or EC number
        if (annotation_product != "") {
            writePfEntry(orf_id, annotation_product, ec_number, start_base, length, output);
        } else if (ec_number != "") {
            // If ec_number valid use original annotation
            annotation_product = fields[9];
            writePfEntry(orf_id, annotation_product, ec_number, start_base, length, output);
        }

        count++;
        if (count % 100 == 0) {
            cout << count << endl;
        }
    }

//    vector <string>::iterator itr;
//    cout << annotations.size() << endl;
//    for(itr = annotations.begin(); itr != annotations.end(); itr++) {
//        cout << *itr << endl;
//    }

    // Print out annotations
//    LIST_NODE *list_itr = my_list->head;
//    while(list_itr != NULL) {
//        cout << list_itr->annotation << endl;
//        list_itr = list_itr->next;
//    }

    input.close();
    output.close();
}

/*
 * Function takes all ingredients for a Pathway Tools pf file entry, creates the string, and writes the string to the
 * .pf output file.
 */
void writePfEntry(string orf_id, string annotation_product, string ec_number, int &start_base, int length, std::ofstream &output) {

    // Example .pf format
    // ID <orf_id>
    // NAME <orf_id>
    // STARTBASE <start_base>
    // ENDBASE <start_base>+<length>
    // PRODUCT-TYPE P
    // FUNCTION <annotation_product>
    // EC <ec_number>
    // METACYC <rxn_id> (optional)
    // PRODUCT-ID <product_id> (optional)
    // //

    string pf_entry = "";
    pf_entry = pf_entry + "ID " + orf_id + "\n";
    pf_entry = pf_entry + "NAME " + orf_id + "\n";
    pf_entry = pf_entry + "STARTBASE " + to_string(start_base) + "\n";
    pf_entry = pf_entry + "ENDBASE " + to_string(start_base + length) + "\n";
    pf_entry = pf_entry + "PRODUCT-TYPE P" + "\n";
    pf_entry = pf_entry + "FUNCTION " + annotation_product + "\n";
    if (ec_number != "")
        pf_entry = pf_entry + "EC " + ec_number + "\n";
    pf_entry = pf_entry + "//" + "\n";
    output << pf_entry;

    // update startbase
    start_base += length + 10;

}

/*
 * Checks an individual annotation (list of words) and its substrings for 'complete' status in in the MetaCyc trie. If
 * complete pointer to position in MetaCyc tree.
 */
string processAnnotationForPtools(vector <char *> annotation_words, PTOOLS_NODE *root, PTOOLS_NODE *ptools_ptr,
                                  bool complete, vector <string> word_list, vector <string> max_word_list, string annotation) {

    // cout << "In processAnnotationForPtools()" << endl;

    // Try to push current word
    for (int i = 0; i < annotation_words.size(); i++) {
        word_list.clear();
        ptools_ptr = root; // reset root
        for (int j = i; j < annotation_words.size(); j++) {
            string word = string(annotation_words[j]);
            if (pushWordForward(word, ptools_ptr)) {
                // Check to see if pointer now at word that completes an annotation
                // cout << "Found " << annotation_words[j] << endl;
                // ptools_ptr.getChildNode(word)->id2;
                word_list.push_back(ptools_ptr->id2);
                if(ptools_ptr->complete) {
                    complete = true;
                    if (word_list.size() > max_word_list.size()) {
                        max_word_list = word_list;
                    }
                }
            } else {
                //cout << "Not Found " << annotation_words[j] << endl;
                break;
            }
        }
    }
    if (complete) {
        // Found complete ptools annotation in annotation_words
        vector <string>::iterator itr;
        for(itr = max_word_list.begin(); itr != max_word_list.end(); itr++) {
            annotation = annotation + *itr + " ";
        }
        return annotation;
        // cout << "Found complete" << endl;
        // cout << "Annotation: " << annotation << endl;
    }
    return "";

}

/*
 * Checks to see if current word is a child of the current node pointed at by ptools_ptr. Moves pointer to child and
 * returns true if so, and returns false otherwise.
 */
bool pushWordForward(string word, PTOOLS_NODE *ptools_ptr) {
    if( ptools_ptr->hasChild(word)) {
        // Word is a child
        ptools_ptr = ptools_ptr->getChildNode(word);
        return true;
    } else {
        return false;
    }
}