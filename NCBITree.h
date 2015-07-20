//
// Created by Niels Hanson on 2015-07-18.
//

#ifndef MP_SPEED_NCBITREE_H
#define MP_SPEED_NCBITREE_H

#include <string>

// #include "types.h"
#include "utilities.h"


using namespace std;

/*
 * Node of the ncbi taxonomy for calculating the Lowest Common Ancestor (LCA) algorithm
 */
struct TREENODE {
    string taxa_id;
    TREENODE* parent;
    bool seen;
    unsigned int count;
    map<string, TREENODE*> children;
    /*
     * Inserts adds a child to the current node
     */
    void insertChild(TREENODE *child) {
        child->parent = this;
        children[child->taxa_id] = child;
    }

    TREENODE() {
        taxa_id = "-1";
        seen = false;
        count = 0;
    }

};

class NCBITree {
    TREENODE *root;
    string ncbi_catalog_file;
    string ncbi_catalog_names_map_file;
    string ncbi_nodes_file;
    public:
        // members
        map<string, string> NCBI_ID_to_Common;
        map<string, string> RefSeqID_to_NCBI_ID;
        map<string, TREENODE*> treeNodeLookup;
        // methods
        NCBITree(string ncbi_catalog, string ncbi_catalog_names_map, string ncbi_nodes); // constructor
        bool BuildCatalogNamesMap();
        bool BuildRefSeqCatalog();
        bool CreateNCBITree();
        string getLCA(vector<string> ncbi_ids);
};

#endif //MP_SPEED_NCBITREE_H
