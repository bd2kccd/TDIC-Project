#include<math.h>
#include <algorithm>
#include <float.h>
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

void Node::calculateCPT(vector<int>& completeMatrix, int numOfCases){
    
    int numOfParents = this->parents.size();
       
    int Aindex = this->index;
    vector<int>& parentIndex = this->parents;

    vector<float> countA1(pow(2,numOfParents),0);
    vector<float> countA0(pow(2,numOfParents),0);
    for(int i = 0; i < pow(2,numOfParents); i++){
        countA1[i] = 0;
        countA0[i] = 0;
    }
    for (int c = 0; c < numOfCases; c++){
        int Avalue = completeMatrix[numOfCases * Aindex + c];

        vector<int> parentValue;
        for (int p = 0; p < numOfParents; p++){
            parentValue.push_back(completeMatrix[numOfCases * parentIndex[p] + c]);
        }

        int countIdx = 0; 
        for (int p = numOfParents - 1; p >= 0; p--){
            countIdx += pow(2,p) * parentValue[p]; 
        }
        if (Avalue == 1)
            countA1[countIdx] ++;
        else
            countA0[countIdx] ++;
    }
    count1.clear();
    count0.clear();
    for (int i = 0; i < countA1.size(); i++){
        count1.push_back(countA1[i]);
        count0.push_back(countA0[i]);
    }
    
    this->CPT1.clear();
    for (int i = 0; i < countA1.size(); i++){
        float c1 = countA1[i];
        float c0 = countA0[i];
        double cpt1 = 0.0;
        if (c1+c0 == 0){
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

void Node::updateACPT(int change, vector<int>& completeMatrix, int numOfCases, int caseNum){
    vector<int>& Aparents = this->parents;
    int numOfParents = Aparents.size();
    
    vector<int> parentValue;
    for (int p = 0; p < numOfParents; p++){
        parentValue.push_back(completeMatrix[numOfCases * Aparents[p] + caseNum]);
    }
    
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
        if (c1+c0 == 0){
            cpt1 = 0.5;
        }else{
            cpt1 = double(c1)/(c1+c0);
        }
        this->CPT1.push_back(cpt1) ;    
    }    
}

void Node::updateChildCPT(int change, int Aindex, vector<int>& completeMatrix, int numOfCases, int caseNum){
    
    vector<int>& cParent = this->parents; 
    int numOfParents = cParent.size();
    if (numOfParents == 0)
        return;
    int cValue = completeMatrix[numOfCases * this->index + caseNum];

    vector<int> parentValuePlus, parentValueMinus;
    for (int p = 0; p < numOfParents; p++){
        int currentValue = completeMatrix[numOfCases * cParent[p] + caseNum];
        parentValuePlus.push_back(currentValue);
        if (cParent[p] == Aindex){
            parentValueMinus.push_back(currentValue - change); 
        }else{
            parentValueMinus.push_back(currentValue);
        }
    }

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
    

    this->CPT1.clear();
    for (int i = 0; i < count1.size(); i++){
        int c1 = count1[i];
        int c0 = count0[i];
        double cpt1 = 0.0;
        if (c1+c0 == 0){
            cpt1 = 0.5;
        }else{
            cpt1 = double(c1)/(c1+c0);
        }
        this->CPT1.push_back(cpt1) ;    
    }    
    
    
}
