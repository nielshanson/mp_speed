#include "idTree.h"
#include <iostream>


using namespace std;

IDTREE::IDTREE() {

    node = new TRIENODE;


}
void IDTREE::createTrie(vector<string> &words) {

    for(vector<string>::iterator it = words.begin(); it != words.end(); it++) {
        this->insert(*it);
    }

}
       
void IDTREE::insert(string word) {
     TRIENODE *node = this->node;
     char c;
     for(unsigned int i = 0; i < word.size(); i++) {
        c = word.at(i); 
        if( node->children.find(c) == node->children.end() ) {
           node->children[c] = new TRIENODE;
           node->children[c]->c = c;
        }
        node= node->children[c];
     }
     node->end = true;
}

string IDTREE::find(string annot) {
    unsigned int a=0, b=1;
    char c;
    TRIENODE *node; 

    while(a < annot.length()) {
         
        if( a> 0 && std::isalnum(annot[a-1])) { a++; continue; }
        
        b = a; 
        node = this->node;
        while( b < annot.length() ) {
           c = annot.at(b); 
           
           if( node->children.find(c) != node->children.end() ) {
              node=node->children[c];
              b++;
              if( node->end && ( b==annot.length() || !std::isalnum(annot[b]) ) ) {
                 return annot.substr(a, b-a);
              }
           }
           else
              break;
        } //inner while

        a++; 
     }
     return string("");
}

void IDTREE::free(TRIENODE *node) {
   
     if(node==0 ) return;

     for(map<char, TRIENODE *>::iterator it = node->children.begin(); it!= node->children.end(); it++ ) {
            this->free(it->second);
     }
     delete node;
}


IDTREE::~IDTREE() {
     TRIENODE *node = this->node; 
     this->free(node);
}

