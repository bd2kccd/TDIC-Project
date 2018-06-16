/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Data.h
 * Author: XIM33
 *
 * Created on April 4, 2018, 2:00 PM
 */
#include <vector>
#include<map>
#include "Node.h"
using namespace std;

#ifndef DATASET_H
#define DATASET_H

class Data {
public:
//    Data();
    Data(string pFileName, string iFileName, string edgeFileName);
    Data(const Data& orig);
    virtual ~Data();
    
    vector<int> &getCombinedMatrix(void) {return combinedMatrix;};
    void inferActivation();
    void calCPTofEachNode();
    double calJointProbOfAllNodes();
    void thresholdNupdate();
    void outputCombinedMatrix(string outPath);
    void outputActivatMatrix(string outPath);
    void outputRandomNumMatrix(string outPath);
    void outputJointProb(string outPath);
private:
    int nProtein; //number of proteins
    int nCase; //number of cases
    
    vector <string> proteinNames;

    vector<int> intervMatrix; //1: stimulator, -1: inhibitor, 0: no intervention
    vector<int> combinedMatrix; //activation nodes + phosphorylation nodes

    vector<double> activatMatrix;
    vector<float> randomNumMatrix;
    
    vector<Node*> nodeList; 
    
    vector<string> edges;
    
    map<string,int> nodeListMap; //nodeName to nodeIndex
    
    vector<double> jointProbs;
    
    
    void readinMatrix(string fileName, char type); //read in phosphorylation table, initial activation table and combine to nodeValueTable 
    void readinEdges(string edgeFileName);    
    void buildNetwork();
    double getCPT1valueOfANode(int Aindex, int caseNumber);
//    double getCPTvalueOfANode(Node* Anode, int caseNumber );
    void getCPTvalueOfOneChildOfANode(int childIndex, int Aindex,int caseNumber, vector<double>& AchildCPT); 
    
    
};

#endif /* DATASET_H */

 