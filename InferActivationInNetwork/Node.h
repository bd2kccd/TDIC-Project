
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
    string getName(void) {return name;};
    void setName(string name ) { this->name = name;};
    
    int getIndex(void) {return index;};
    void setIndex(int idx){this->index = idx;};
    
    char getType(void) {return type;};
    void setType(char type){this->type = type;};
    
    vector<int>& getParents(void) {return this->parents;};
    void addToParents(int nodeIndex);
    
    vector<int>& getChildren(void) {return this->children;};
    void addToChildren(int nodeIndex);
    
    vector<double>& getCPT(){return this->CPT1;};
    void calculateCPT(vector<int>& combinedMatrix, int numOfCases);
    void updateACPT(int change, vector<int>& combinedMatrix, int numOfCases, int caseNum);
    void updateChildCPT(int change,int parentIdx, vector<int>& combinedMatrix, int numOfCases, int caseNum);
    
private:
    string name;
    int index;
    char type;
    vector<int> parents; 
    vector<int> children;
    vector<double> CPT1; 
    vector<float> count1, count0; 
};

#endif /* NODE_H */

