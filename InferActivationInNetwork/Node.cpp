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
#include <float.h>
//#include<iostream>
#include "Node.h"

using namespace std;

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

void Node::calculateCPT(vector<int>& combinedMatrix, int numOfCases){
    
    int numOfParents = this->parents.size();
       
    int Aindex = this->index;
    vector<int>& parentIndex = this->parents;
    //loop through cominedMatrix to count parent nodes combination for A=1 and A=0
    vector<float> countA1(pow(2,numOfParents),0);
    vector<float> countA0(pow(2,numOfParents),0);
    for(int i = 0; i < pow(2,numOfParents); i++){
        countA1[i] = 0.1;
        countA0[i] = 0.1;
    }
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
    count1 = countA1;
    count0 = countA0;
    //calculate CPT1 (we only save p(A=1|.), P(A=0|.)=1-p(A=1|.)
    this->CPT1.clear();
    for (int i = 0; i < countA1.size(); i++){
        float c1 = countA1[i];
        float c0 = countA0[i];
        double cpt1 = 0.0;
        if (c1+c0 == 0){//the combination never appear in real data
            cpt1 = 0.5;
        }else{
            cpt1 = double(c1)/(c1+c0);
        }
        this->CPT1.push_back(cpt1) ;    
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

void Node::updateACPT(int change, vector<int>& combinedMatrix, int numOfCases, int caseNum){
    vector<int>& Aparents = this->parents;
    int numOfParents = Aparents.size();
    
    //get value of A's parents
    vector<int> parentValue;
    for (int p = 0; p < numOfParents; p++){
        parentValue.push_back(combinedMatrix[numOfCases * Aparents[p] + caseNum]);
    }
    //calculate the count index
    int countIdx=0; 
    for (int p = numOfParents - 1; p >= 0; p--){
        countIdx += pow(2,p) * parentValue[p];
    }
    if (change == 1){ //0->1
        count1[countIdx] ++;
        count0[countIdx] --;
    }
    else{//1->0
        count0[countIdx] ++;
        count1[countIdx] --;
    }

    this->CPT1.clear();
    for (int i = 0; i < count1.size(); i++){
        float c1 = count1[i];
        float c0 = count0[i];
        double cpt1 = 0.0;
        if (c1+c0 == 0){//the combination never appear in real data
            cpt1 = 0.5;
        }else{
            cpt1 = double(c1)/(c1+c0);
        }
        this->CPT1.push_back(cpt1) ;    
    }    
}

void Node::updateChildCPT(int change, int Aindex, vector<int>& combinedMatrix, int numOfCases, int caseNum){
    //change = 1: 0->1; change = -1: 1->0
    vector<int>& cParent = this->parents; 
    int numOfParents = cParent.size();
    if (numOfParents == 0)
        return;
    //for test purpose
    int pos = -1;
    for (int i = 0; i < numOfParents; i++){
        if (cParent[i] == Aindex){
            pos = i;
            break;
        }
    }
    if(pos == -1){
//        cout << "error: parent index is not found in undateChildCPT.";
        exit(1);
    }
    //test end
    int cValue = combinedMatrix[numOfCases * this->index + caseNum];
    //get value of parents of A node
    vector<int> parentValuePlus, parentValueMinus;
    for (int p = 0; p < numOfParents; p++){
        int currentValue = combinedMatrix[numOfCases * cParent[p] + caseNum];
        parentValuePlus.push_back(currentValue);
        if (cParent[p] == Aindex){
            parentValueMinus.push_back(currentValue - change); //if 0 the 1, if 1 then 0
        }else{
            parentValueMinus.push_back(currentValue);
        }
    }
    //calculate the count index
    int countIdxPlus = 0,countIdxMinus = 0; 
    for (int p = numOfParents - 1; p >= 0; p--){
        countIdxPlus += pow(2,p) * parentValuePlus[p];
        countIdxMinus += pow(2,p) * parentValueMinus[p];
    }
    if (cValue == 1){
        count1[countIdxPlus] ++;
        count1[countIdxMinus] --;
    }
    else{
        count0[countIdxPlus] ++;
        count0[countIdxMinus] --;
    }
    

    //calculate CPT1 (we only save p(A=1|.), P(A=0|.)=1-p(A=1|.)
    this->CPT1.clear();
    for (int i = 0; i < count1.size(); i++){
        int c1 = count1[i];
        int c0 = count0[i];
        double cpt1 = 0.0;
        if (c1+c0 == 0){//the combination never appear in real data
            cpt1 = 0.5;
        }else{
            cpt1 = double(c1)/(c1+c0);
        }
        this->CPT1.push_back(cpt1) ;    
    }    
    
    
}


//void Node::sortParentsChildren(){
//    std::sort(this->parents.begin(), this->parents.end());
//    std::sort(this->children.begin(), this->children.end());
//}