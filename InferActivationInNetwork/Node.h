/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Node.h
 * Author: XIM33
 *
 * Created on April 9, 2018, 3:18 PM
 */

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
    void calculateCPT(vector<unsigned int>& combinedMatrix, int numOfNodes, int numOfCases);
    
    void sortParentsChildren();
    
private:
    string name;
    int index;
    char type;
    vector<int> parents; //save parent node index
    vector<int> children;//using Node itself to easy get children's parents
    vector<double> CPT1; //P(N-1|Pa(N))  for   P(N-0|Pa(N)) = 1-P(N-1|Pa(N))
};

#endif /* NODE_H */

