//
// Created by Niels Hanson on 2015-07-18.
//

#ifndef MP_SPEED_NCBITREE_H
#define MP_SPEED_NCBITREE_H

#include <string>

#include "utilities.h"
#include "types.h"

using namespace std;

class NCBITree {
    TREENODE *root;
    string ncbi_catalog_file;
    string ncbi_catalog_names_map_file;
    string ncbi_nodes_file;
    map<string, string> NCBI_ID_to_Common;
    map<string, string> RefSeqID_to_NCBI_ID;
    map<string, TREENODE*> treeNodeLookup;
    public:
        NCBITree(string ncbi_catalog, string ncbi_catalog_names_map, string ncbi_nodes); // constructor
        bool BuildCatalogNamesMap();
        bool BuildRefSeqCatalog();
        bool CreateNCBITree();
        string getLCA(vector<string> ncbi_ids);
};

#endif //MP_SPEED_NCBITREE_H
