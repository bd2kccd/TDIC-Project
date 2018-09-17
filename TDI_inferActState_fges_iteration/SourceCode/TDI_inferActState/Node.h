
#include <vector>
#include <set>
#include <string>
using namespace std;

#ifndef NODE_H
#define NODE_H

class Node {
public:
    Node();
    Node(string nodeName, char nodeType, int nodeIndex);
    Node(const Node& orig);
    virtual ~Node();
    string getName(void) {return this->name;};
    void setName(string name ) { this->name = name;};
    
    int getIndex(void) {return this->index;};
    void setIndex(int idx){this->index = idx;};
    
    char getType(void) {return this->type;};
    void setType(char type){this->type = type;};
    
    vector<int>& getParents(void) {return this->parents;};
    void addToParents(int nodeIndex);
    
    vector<int>& getChildren(void) {return this->children;};
    void addToChildren(int nodeIndex);
    
    vector<double>& getCPT(){return this->CPT1;};
    void calculateCPT(vector<int>& completeMatrix, int numOfCases);
    void updateACPT(int change, vector<int>& completeMatrix, int numOfCases, int caseNum);
    void updateChildCPT(int change,int parentIdx, vector<int>& completeMatrix, int numOfCases, int caseNum);
    
private:
    string name;
    int index;
    char type;
    vector<int> parents; 
    vector<int> children;
    vector<double> CPT1; //we only save p(A=1|.), P(A=0|.)=1-p(A=1|.)
    vector<float> count1, count0; 
};

#endif /* NODE_H */

