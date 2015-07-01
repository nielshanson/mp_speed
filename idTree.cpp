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
     for(int i = 0; i < word.size(); i++) {
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
     unsigned int a=0, b=1, d;
     char c;
     TRIENODE *node; 

     while(a < annot.length()) {
        b = a; 
        node = this->node;
        while( b < annot.length() ) {
           c = annot.at(b); 

           if( node->children.find(c) != node->children.end() ) {
              node=node->children[c];
              b++;
              if( node->end && ( b==annot.length() || !std::isalnum(annot[b]) ) ) {
                 return annot.substr(a, b-a+1);
              }
           }
           else 
              break;
        } //inner while

        a++; 
     }
     return string("");
}
