/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Node.cpp
 * Author: XIM33
 * 
 * Created on April 9, 2018, 3:18 PM
 */
#include<math.h>
#include <algorithm>
#include "Node.h"

Node::Node() {
}

Node::Node(string nodeName, char nodeType, int nodeIndex) {
    this->name = nodeName;
    this->type = nodeType;
    this->index = nodeIndex;
}

Node::Node(const Node& orig) {
}

Node::~Node() {
}

void Node::calculateCPT(vector<unsigned int>& combinedMatrix, int numOfNodes, int numOfCases){
    
    int numOfParents = this->parents.size();
       
    int Aindex = this->index;
    vector<int>& parentIndex = this->parents;
    //loop through cominedMatrix to count parent nodes combination for A=1 and A=0
    //first need to calculate the count index = 2^(numOfParent-1)*P1value + 2^(numOfParent-2)*P2value + ... + 2^0 * P(numOfParent)value(last one)
    vector<int> countA1(pow(2,numOfParents),0);
    vector<int> countA0(pow(2,numOfParents),0);
    for (int c = 0; c < numOfCases; c++){
        int Avalue = combinedMatrix[numOfCases * Aindex + c];
        //get value of parents of A node
        vector<int> parentValue;
        for (int p = 0; p < numOfParents; p++){
            parentValue.push_back(combinedMatrix[numOfCases * parentIndex[p] + c]);
        }
        //calculate the count index
        int countIdx = 0; 
        for (int p = numOfParents - 1; p >= 0; p--){
            countIdx += pow(2,p) * parentValue[p]; 
        }
        if (Avalue == 1)
            countA1[countIdx] ++;
        else
            countA0[countIdx] ++;
    }
    
    //calculate CPT1 (we only save p(A=1|.), P(A=0|.)=1-p(A=1|.)
    this->CPT1.clear();
    for (int i = 0; i < countA1.size(); i++){
        int c1 = countA1[i];
        int c0 = countA0[i];
        this->CPT1.push_back(double(c1)/(c1+c0)) ;           
    }
}

void Node::addToParents(int nodeIndex){
    if(std::find(this->parents.begin(), this->parents.end(), nodeIndex) == this->parents.end()) {
        this->parents.push_back(nodeIndex);
    }
}

void Node::addToChildren(int nodeIndex){
    if(std::find(this->children.begin(), this->children.end(), nodeIndex) == this->children.end()) {
        this->children.push_back(nodeIndex);
    }
}

void Node::sortParentsChildren(){
    std::sort(this->parents.begin(), this->parents.end());
    std::sort(this->children.begin(), this->children.end());
}