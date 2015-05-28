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

    cout << options.ptools_rxn_file << endl;

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
            if(! start->hasChild(words[i])) {
                // create child node
                start->insertChild(words[i]);
                if (i == (words.size()-1)) {
                    // flag node as finished
                    start->children[words[i]]->complete = true;
                }
            }
            start = start->children[words[i]];
        }
    }
    cout << "Loaded tree" << endl;

    // print out words in ptools tree
    start = root;
    vector<PTOOLS_NODE *> node_stack;
    vector<string> pstack;
    vector<string> path;

    node_stack.push_back(start);
    string headpath = "";
    string path_line = "";
    path.push_back(headpath);
    while(!node_stack.empty()) {
        PTOOLS_NODE* top = node_stack.back();
        node_stack.pop_back();
        headpath = path.back(); // keep track of current path
        path.pop_back();
        if (top->children.size() == 0) {
            // print the stack
            // cout << "headpath: " << headpath << endl;
            // cout << "pstack: " << pstack.back() << endl;
            path_line = headpath;
            cout << path_line << endl;

        }
        map<string, PTOOLS_NODE*> children = top->children;
        for (std::map<string, PTOOLS_NODE*>::iterator itr=children.begin(); itr != children.end(); ++itr) {
            node_stack.push_back(itr->second);
            path.push_back(headpath + "->"+ itr->first);
        }

    }
    exit(1);

    string test_anno = "The cat in the hat";

    // Test out linked list
    LIST *my_list = new LIST(root);

    processAnnotationsForPTools(my_list, root, options.annotation_table);

    exit(1);
    // insert into the list
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


    exit(1);

}

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

    // Get header
    while( std::getline(input, line ).good()) {
        // skip header
        if (count == 0) {
            count++;
            continue;
        }

        // std::cout << line << std::endl;
        split(line, fields, buf,'\t');
        ptools_list.push_back(fields[1]);
        //if( line.size()==0 or  line[0]!='>') continue;

        //split_seq_name(line, fields, this->buf);

//        name = (fields[0]);
//
//        //  std::cout << name << std::endl;
//        //std::cout << fields[1] << std::endl;
//
//        if( fields.size()< 2)
//            annotation = "hypothetical protein";
//        else
//            annotation = string(fields[1]);
//
//        if( query_dictionary.find(name) != query_dictionary.end() )
//            annot_map->insert(std::make_pair(name,annotation));
//
//        if (count%PRINT_INTERVAL==0)
//            std::cout << count << std::endl;
        count++;
    }

    // std::cout << "Number of annotatons scanned " << count << std::endl;
    // std::cout << "Number of annotation loaded " <<  annot_map->size() << std::endl;

    input.close();
}

void processAnnotationsForPTools(LIST *my_list, PTOOLS_NODE *root, string annotation_file) {
    // For each annotation in annotation_file
    std::cout << "Reading annotation_file " << annotation_file <<  std::endl;
    input.open(annotation_file.c_str(), std::ifstream::in);
    if(!input.good()){
        std::cerr << "Error opening '"<< annotation_file <<"'. Bailing out." << std::endl;
        return ;
    }

    // Scanning variables
    int count =0;
    string line;
    char buf[10000]; // Temp buffer
    vector <char *> fields; // Vector for parsed fields
    vector <char *> annotation_words; // Vector of annotation words

    while( std::getline(input, line ).good()) {

        // cout << line << endl;
        split(line, fields, buf, '\t');
        split(fields[9], annotation_words, buf, ' ');

        processAnnotationForPtools(annotation_words, root, my_list);

        //if( line.size()==0 or  line[0]!='>') continue;

        //split_seq_name(line, fields, this->buf);

//        name = (fields[0]);
//
//        //  std::cout << name << std::endl;
//        //std::cout << fields[1] << std::endl;
//
//        if( fields.size()< 2)
//            annotation = "hypothetical protein";
//        else
//            annotation = string(fields[1]);
//
//        if( query_dictionary.find(name) != query_dictionary.end() )
//            annot_map->insert(std::make_pair(name,annotation));
//
//        if (count%PRINT_INTERVAL==0)
//            std::cout << count << std::endl;
        count++;
    }

    // std::cout << "Number of annotatons scanned " << count << std::endl;
    // std::cout << "Number of annotation loaded " <<  annot_map->size() << std::endl;

    input.close();

}

void processAnnotationForPtools(vector <char *> annotation_words, PTOOLS_NODE *root, LIST *my_list) {
    PTOOLS_NODE ptools_ptr = *root;
    bool complete = false;
    // Try to push current word
    for (int i = 0; i < annotation_words.size(); i++) {
        for (int j = i; j < annotation_words.size(); j++) {
            if (pushWordForward(annotation_words[j], ptools_ptr)) {
                // check to see if pointer now at word that completes an annotation
                // cout << ptools_ptr->id1 << endl;
                cout << "Found " << annotation_words[j] << endl;
                if(ptools_ptr.complete) {
                    complete = true;
                }
            } else {
                cout << "Not Found " << annotation_words[j] << endl;
                break;
            }
        }
    }
    if (complete) {
        // found complete ptools annotation in annotation_words
        cout << "Found complete" << endl;
        my_list->insert(&ptools_ptr);
    }

}

bool pushWordForward(char *word, PTOOLS_NODE &ptools_ptr) {
    if( ptools_ptr.children.find(word) != ptools_ptr.children.end()) {
        // Word is a child
        ptools_ptr = *ptools_ptr.children[word];
        return true;
    } else {
        return false;
    }


}