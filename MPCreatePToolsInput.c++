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
            }
            start = start->children[words[i]];
        }
    }
    cout << "Loaded tree" << endl;

    string test_anno = "The cat in the hat";

    // Test out linked list
    LIST *my_list = new LIST(root);

    // processAnnotationsForPTools(my_list, ptools_tree);

    // insert into the list
    my_list->insert("One",1);
    my_list->insert("Two",1);
    my_list->insert("Three",1);

    // Move current to next node
    my_list->nextNode();
    // Test delete
    my_list->deleteCurr();
    // Test insert
    my_list->insertAtCurr("Four",1);

    // Print out list
    LIST_NODE *list_itr = my_list->head;
    while(list_itr != NULL) {
        cout << list_itr->my_data << endl;
        list_itr = list_itr->next;
    }


    // print out ptools tree
    // start = head;
    // print_dfs(start, "");

    // Create listnode linked list


    exit(1);

    // print out words in ptools tree
//    start = head;
//    stack<PTOOLS_NODE *> node_stack;
//    stack<string> pstack;
//
//    node_stack.push(start);
//    pstack.push("");
//    string line = "";
//    while(!node_stack.empty()) {
//        pstack.pop();
//        PTOOLS_NODE* top = node_stack.top();
//        node_stack.pop();
//        if (top->children.size() == 0) {
//            line = "";
//            while(!pstack.empty()) {
//                line = line + " " + pstack.top();
//                pstack.pop();
//            }
//            cout << line << endl;
//            line = "";
//        }
//        map<string, PTOOLS_NODE*> children = top->children;
//        cout << "Children: " << children.size() << endl;
//        for (std::map<string, PTOOLS_NODE*>::iterator itr=children.begin(); itr != children.end(); ++itr) {
//            pstack.push(itr->first);
//            node_stack.push(itr->second);
//            cout << itr->first << endl;
//        }
//    }

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