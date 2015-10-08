//
// Created by Niels Hanson on 2015-07-18.
//

#include "NCBITree.h"

NCBITree::NCBITree(string ncbi_catalog_file = "",
                   string ncbi_catalog_names_map_file = "",
                   string ncbi_nodes_file = "",
                   string megan_map_file = "") {
    this->ncbi_catalog_file = ncbi_catalog_file;
    this->ncbi_catalog_names_map_file = ncbi_catalog_names_map_file;
    this->ncbi_nodes_file = ncbi_nodes_file;
    this->megan_map_file = megan_map_file;
    this->built = false;
}

bool NCBITree::init() {
    this->BuildRefSeqCatalog();
    this->BuildCatalogNamesMap();
    this->CreateNCBITree();
    this->built = true;
    return true;
}

bool NCBITree::BuildRefSeqCatalog() {
    ifstream input;
    string filename = this->ncbi_catalog_names_map_file;
    string line;
    char buf[1000]; // Temp buffer
    vector<char *> fields; // Vector for parsed fields

    cout << "Building NCBI to Common Name Map:" << endl;

    input.open(filename.c_str(), std::ifstream::in);
    if(!input.good()){
        cerr << "Error opening '"<< filename << "'. " << endl;
        return false;
    }

    this->NCBI_ID_to_Common["-1"] = "unknown taxonomy"; // no match

    // split NCBI_ID and common_name
    while( std::getline( input, line ).good() ) {
        split(line, fields, buf, '\t');
        if(this->NCBI_ID_to_Common.find(fields[0]) == this->NCBI_ID_to_Common.end()) {
            this->NCBI_ID_to_Common[fields[0]] = fields[1];
        }
    }
    input.close();

    if (this->megan_map_file != "") {
        filename = this->megan_map_file;
        fields.clear();
        input.open(filename.c_str(), std::ifstream::in);
        if(!input.good()){
            cerr << "Error opening '"<< filename << "'. " << endl;
            return false;
        }
        // split NCBI_ID and common_name
        while( std::getline( input, line ).good() ) {
            split(line, fields, buf, '\t');
            if(this->NCBI_ID_to_Common.find(fields[0]) == this->NCBI_ID_to_Common.end()) {
                this->NCBI_ID_to_Common[fields[0]] = fields[1];
            }
        }
        input.close();
    }

    return true;
}

/*
 * Takes ncbi_catalog file and splits fields
 */
bool NCBITree::BuildCatalogNamesMap() {
    ifstream input;
    string filename = this->ncbi_catalog_file;
    string line;
    char buf[1000]; // Temp buffer
    vector<char *> fields; // Vector for parsed fields

    input.open(filename.c_str(), std::ifstream::in);
    if(!input.good()){
        cerr << "Error opening '"<< filename << "'. " << endl;
        return false;
    }

    int count = 0;
    // split RefSeq_ID and NCBI_ID
    while( std::getline( input, line ).good() ) {
        split(line, fields, buf, '\t');
        // NCBI_ID second field
        count++;
        this->RefSeqID_to_NCBI_ID[fields[1]] = fields[0];
        if (count % 1000000 == 0) cout << "Reading RefSeqID to NCBI_ID: " << count << endl;
    }
    input.close();

    return true;
}

bool NCBITree::CreateNCBITree() {
    ifstream input;
    string filename = this->ncbi_nodes_file;
    string line;
    char buf[1000]; // Temp buffer
    vector<char *> fields; // Vector for parsed fields

    input.open(filename.c_str(), std::ifstream::in);
    if(!input.good()){
        cerr << "Error opening '"<< filename << "'. " << endl;
        return false;
    }

    int count = 0;

    // Declare memory for child and parent nodes
    TREENODE *child_node;
    TREENODE *parent_node;

    // split RefSeq_ID and NCBI_ID
    while( std::getline( input, line ).good() ) {
        split(line, fields, buf, '\t');
        // NCBI_ID second field
        count++;
        string child = fields[1];
        string parent = fields[2];

        if (child != parent) {
            // Create child node (if needed)
            if (this->treeNodeLookup.find(child) == this->treeNodeLookup.end()) {
                child_node = new TREENODE();
                child_node->taxa_id = child;
                this->treeNodeLookup[child] = child_node;
            } else {
                child_node = treeNodeLookup[child];
            }
            // Create parent node (if needed)
            if (this->treeNodeLookup.find(parent) == this->treeNodeLookup.end()) {
                parent_node = new TREENODE();
                parent_node->taxa_id = parent;
                this->treeNodeLookup[parent] = parent_node;
            } else {
                parent_node = treeNodeLookup[parent];
            }

            parent_node->insertChild(child_node);
        }
        if (count % 1000000 == 0) cout << "Creating NCBI Taxononmy Tree: " << count << endl;
    }
    input.close();

    // set root
    this->root = treeNodeLookup["1"];

    return true;
}

vector<string> NCBITree::getLineage(string ncbi_id) {
    TREENODE *cur_node;
    vector<string> lineage_ids = vector<string>();
    if (this->treeNodeLookup.find(ncbi_id) != this->treeNodeLookup.end()) {
        cur_node = this->treeNodeLookup[ncbi_id];
        while(cur_node != 0) {
            lineage_ids.push_back(cur_node->taxa_id);
            cur_node = cur_node->parent;
        }
    }
    return lineage_ids;
}

string NCBITree::getLCA(vector<string> ncbi_ids) {
    vector<string> ncbi_ids_clean; // confirmed ids
    for(vector<string>::iterator itr = ncbi_ids.begin(); itr != ncbi_ids.end(); ++itr) {
        if (this->treeNodeLookup.find(*itr) != this->treeNodeLookup.end()) {
            ncbi_ids_clean.push_back(*itr);
        }
    }

    string lca = "-1"; // no answer

    if (ncbi_ids_clean.size() == 0) {
        return lca;
    } else if (ncbi_ids_clean.size() == 1) {
        return ncbi_ids_clean[0];
    } else {
        bool lca_flag;
        TREENODE *cur_node;
        for(vector<string>::iterator itr = ncbi_ids_clean.begin(); itr != ncbi_ids_clean.end(); ++itr) {
            cur_node = this->treeNodeLookup[*itr];
            lca_flag = false;
            while (cur_node != 0) {
                cur_node->count++;
                if (!cur_node->seen) {
                    cur_node->seen = true;
                } else {
                    if (!lca_flag && cur_node->count == ncbi_ids_clean.size()) {
                        lca = cur_node->taxa_id;
                        lca_flag = true;
                    }
                }
                cur_node = cur_node->parent;
            }
        }

        // reset lca tree
        for(vector<string>::iterator itr = ncbi_ids_clean.begin(); itr != ncbi_ids_clean.end(); ++itr) {
            cur_node = this->treeNodeLookup[*itr];
            while (cur_node != 0) {
                cur_node->seen = false;
                cur_node->count = 0;

                cur_node = cur_node->parent;
            }

        }
    }

    return lca;
}

