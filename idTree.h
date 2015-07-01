#ifndef _IDTREE
#define _IDTREE

#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;


typedef struct _TRIENODE TRIENODE;

typedef struct _TRIENODE  {
     _TRIENODE() {  end = false; }

     char c;
     map<char, TRIENODE *> children;
     bool end;

} TRIENODE;

class IDTREE {

    public:
        IDTREE();
        ~IDTREE();
        void createTrie(vector<string> &words);
        string find(string word);

    private:
        void free(TRIENODE *node);
        TRIENODE *node;
        void insert(string word);

};
#endif
