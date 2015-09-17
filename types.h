#ifndef __PARSE_TYPE
#define __PARSE_TYPE
#include <utility>
#include <vector>
#include <iostream>
#include <map>
#include <stdio.h>

#include "MPParseBlastOptions.h"
#include "MPAnnotateOptions.h"
#include "utilities.h"
#include "idTree.h"
#include "NCBITree.h"

typedef struct _GLOBAL_PARAMS {
    float lambda;
    float k;
} GLOBAL_PARAMS;

typedef struct _BLASTOUT_DATA BLASTOUT_DATA;

typedef struct _THREAD_DATA {
    vector<string>  lines;
    vector<BLASTOUT_DATA>  output[2];
    vector< std::pair<std::string, std::string> > refscorePairs;
    map<std::string, float> refscores;
    float ln2, lnk, lambda;
    MPParseBlastOptions options;
    map<std::string, std::string> *annot_map;
    unsigned int b;
} THREAD_DATA;

/*
 * Structure to assist in the parsing of *.tree.txt files
 */ 
struct HNODE {
    string name;
    string alias;
    int depth;
    vector<HNODE*> children;
};

typedef enum {KEGG, COG, SEED, CAZY} DBTYPE;

/*
 * Structure to keep track of databases, file names, database weights, database types ( KEGG, COG, SEED, etc.)
 * and create temporary files for writing out the functional_and_taxonomic_table.txt contining the annotations
 * for each ORF
 */
typedef struct _DB_INFO {
    vector<string> db_names;
    vector<string> db_tags;
    vector<string> input_blastouts;
    vector<float> weight_dbs;
    vector<DBTYPE> dbtypes;
    map<string, string (*)(const char *) > idextractors;

    string TempFile1(const string &sample, const string &db)  {
        string file_name;
        file_name = sample + string(".") + db + string(".") + string("DB_annotation_table.txt.tmp1");
        return file_name;
    }

    string TempFile2(const string &sample, const string &db)  {
        string file_name;
        file_name = sample + string(".") + db + string(".") + string("DB_annotation_table.txt.tmp2");
        return file_name;
    }

    string FinalFile(const string &sample, const string &db)  {
        string file_name;
        file_name = sample + string(".") + db + string(".") + string("DB_annotation_table.txt");
        return file_name;
    }

    void RemoveTempFile1(const string &sample, const string &db)  {
        string file_name;
        file_name = sample + string(".") + db + string(".") + string("DB_annotation_table.txt.tmp1");
        remove(file_name.c_str());
    }

    void  RemoveTempFile2(const string &sample, const string &db)  {
        string file_name;
        file_name = sample + string(".") + db + string(".") + string("DB_annotation_table.txt.tmp2");
        remove(file_name.c_str());
    }

} DB_INFO;

typedef struct _BLASTOUT_DATA {
    string query, target, product, ec, taxonomy;
    unsigned int q_length, aln_length ;
    float bitscore, bsr, expect, identity;
    void setDefault() {
        query = "";
        target = "";
        product = "";
        ec="";

        q_length = 0;
        aln_length = 0;
        bitscore = 0.0;
        bsr = 0.0;
        expect = 100.0;
        identity = 0.0;
    }
} BLASTOUT_DATA;

typedef struct _WRITER_DATA {
    THREAD_DATA *thread_data;
    unsigned int num_threads;
    std::ofstream output;
} WRITER_DATA;


/* Annotation data structure containing all the data from a parsed B/LASTout 
 * annotation line.
 */
typedef struct _ANNOTATION {
    // functional_and_taxonomic fields: ORF_ID	ORF_length	start	end	Contig_Name	Contig_length	strand	ec	taxonomy	product
    float bsr, value, bitscore;
    string ec, product, taxonomy;
    string query, target;
    
    string orf_id, contig_name, strand;
    unsigned int start, end, length;
    map<string, string> db_ids;
    bool ptools_match; // flag to say annotation comes from trie
    string dbname; // database origin
    unsigned int annotation_score;
} ANNOTATION;


typedef map<string, map<string, vector<ANNOTATION *> > > ANNOTATION_RESULTS;

typedef vector<string> DB_HIT;

/*
 * PTOOLS_NODE that makes up the pathway tools annotation tree
 */
struct PTOOLS_NODE {

    string id1; // lowercase
    string id2; // fullcase
    map<string, PTOOLS_NODE*> children;
    // int count; // annotation count
    bool complete; // Marks that this node ends a complete strand of annotations

    /*
     * Checks to see if this node has a child with a particular name w
     */
    bool hasChild(string w) {
        string lower =  to_lower(w);
        if (children.find( lower ) != children.end()) {
            return true;
        }
        return false;
    }

    /*
     * Inserts adds a child to the current node
     */
    void insertChild(string w) {
        PTOOLS_NODE *child = new PTOOLS_NODE; // create child
        string lower = to_lower(w);
        child->id1 = lower;
        child->id2 = w;
        child->complete = false; // default not leaf
        children[lower] = child;
    }
    /*
     * Sets current node to finished.
     */
    void flagChildFinished(string w) {
        children[to_lower(w)]->complete = true;
    }
    /*
     * Gets the child node of a particular annotation string.
     */
    PTOOLS_NODE* getChildNode(string w) {
        string lower = to_lower(w);
        if (this->hasChild(lower)) {
            return children[lower];
        } else {
            return this;
        }

    }
};


// Data structure for MPAnnotateThreads
typedef struct _THREAD_DATA_ANNOT {
    vector<ANNOTATION>  lines; // annotation lines
    ANNOTATION_RESULTS annot_objects;
    vector<ANNOTATION>  output[2]; 
    vector<string> orfids; // orf_ids this particular thread deals with
    vector<short int> annot_from_db;
    map<string, ANNOTATION> db_hits;
    map<string, IDTREE*> dbNamesToHierarchyIdTree;
    // map<string, map<string, string> > dbNamesToHierarchyIdentifierMaps;
    map<string, map<string, int> > dbNamesToHierachyIdentifierCounts;
    vector<ANNOTATION> metaCycHits;
    PTOOLS_NODE *root;
    NCBITree *ncbi_tree; // formatted NCBI hierarchy for calculating LCA

    MPAnnotateOptions options;
    unsigned int b;
    unsigned int tot_annot;
    map<string, unsigned int> dbwise_annot_count;
    DB_INFO db_info;

    _THREAD_DATA_ANNOT() {
        tot_annot=0;
    }

    void clear() {
        lines.clear();
        dbNamesToHierachyIdentifierCounts.clear();
        metaCycHits.clear();
        // annot_object.clear();
        orfids.clear();
        // annot_from_db.clear();
        // db_hits.clear();
    }

} THREAD_DATA_ANNOT;


// Data structure for writer thread
typedef struct _WRITER_DATA_ANNOT {
    THREAD_DATA_ANNOT  *thread_data;
    unsigned int num_threads;
    DB_INFO db_info;
    map<string, map<string, int> > globalDbNamesToHierachyIdentifierCounts;
    map<string, ANNOTATION > globalMetaCycNamesToAnnotations;
    map<string, map<string, int> > globalMetaCycNamesToDbCounts;
    MPAnnotateOptions options;
    
} WRITER_DATA_ANNOT;


// MPCreatePToolsInput functions

/*
 * LIST_NODE forms an linked list of pointers into nodes of the Pathway Tools annotation tree
 */
struct LIST_NODE {
    LIST_NODE *next, *prev; // pointers to the next and previous nodes
    int index; // index
    string my_data; // some data to hold
    string annotation; // pointer into pathway_tools annotation tree

    LIST_NODE(string annotation) {
        this->annotation = annotation;
    }
};

/*
 * Linked list of LIST_NODEs that point to PTOOLS_NODES in the PTOOLS annotation tree.
 */
struct LIST {
    LIST_NODE *curr, *head; // Current pointer, head pointer
    PTOOLS_NODE *ptools_root; // root of the PTOOLS_TREE

    // Constructor
    LIST(PTOOLS_NODE *my_ptools_root){
        ptools_root = my_ptools_root; // set reference to root of PTOOLS_TREE
        curr = NULL;
        head = NULL;
    }

    // Inserts a node at the head of the linked list
    void insert(string annotation) {
        LIST_NODE *my_node = new LIST_NODE(annotation); // Create Node
        if (head == NULL) {
            head = my_node;
            curr = my_node;
        } else {
            // insert node to head of list
            head->prev = my_node;
            my_node->next = head;
            head = my_node;
            curr = my_node;
        }
    };
    // Deletes the node pointed to by the head pointer
    void remove() {
        if (head == NULL) {
            // Nothing to remove
            return;
        }
        LIST_NODE *del = head;
        head = head->next;
        delete del;
    };
    // Move the current pointer to the next node
    void nextNode() {
        if (curr == NULL) {
            return;
        }
        else {
            if (curr->next != NULL) {
                curr = curr->next;
            }
        }
    };
    // Move the current pointer to the previous node
    void prevNode() {
        if (curr == NULL) {
            return;
        }
        else {
            if (curr->prev != NULL) {
                curr = curr->prev;
            }
        }
    };
    // Delete the node at the current pointer
    void deleteCurr(){
        if (curr == NULL) {
            // Nothing to remove
            return;
        }
        LIST_NODE *del = curr; // node to delete
        if (curr->next == NULL &&  curr->prev != NULL) {
            // Current pointer at end
            curr->prev->next = NULL;
            delete del;
        }
        if (curr->prev == NULL && curr->next != NULL) {
            // Current pointer at beginning
            curr->next->prev = curr->prev; // i.e., NULL
            head = curr->next; // make sure to update head
            curr = curr->next;
            delete del;
        }
        if (curr->next != NULL && curr->prev != NULL) {
            // Current pointer in middle
            curr->prev->next = curr->next;
            curr->next->prev = curr->prev;
            curr = curr->prev;
            delete del;
        }
    }
    // Inserts a node at the current pointer
    void insertAtCurr(string annotation) {
        LIST_NODE *my_node = new LIST_NODE(annotation); // Create Node
        if (curr == NULL) {
            head = my_node;
            curr = my_node;
        } else {
            // insert node after node currently being pointed to by current pointer
            my_node->prev = curr;
            my_node->next = curr->next;
            curr->next = my_node;
            // move current pointer to my_node just inserted
            curr = my_node;
        }
    }
};


#endif



