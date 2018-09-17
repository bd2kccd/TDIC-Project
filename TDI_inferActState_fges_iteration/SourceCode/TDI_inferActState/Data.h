#include <vector>
#include<map>
#include "Node.h"
using namespace std;

#ifndef DATASET_H
#define DATASET_H


class Data {
public:

    Data(string MatrixFileName, string SGAFileName, string DEGFileName, string edgeFileName);
//    Data(const Data& orig);
    virtual ~Data();
    
    vector<int> &getCombinedMatrix(void) {return combinedMatrix;};
    void getInferState();
    void calCPTofEachNode();
    double calJointProbOfAllNodes();
    void outputCombinedMatrix(string outPath);
    void outputActivatMatrix(string outPath);
    void outputJointProb(string outPath);
    void outputCPT(string outpath);
    void convertToCombination(vector<int> &combin_arr, unsigned int n);
private:
    int nDriver; //number of SGA
    int nNode; //number of SGA+DEG
    int nTumor; //number of 
    
    vector <string> nodeNames;

    vector<int> intervMatrix; //1: stimulator, -1: inhibitor, 0: no intervention
    vector<int> combinedMatrix; //activation nodes + phosphorylation nodes

    vector<double> activatMatrix;
    
    vector<Node*> nodeList; 
    
    vector<string> edges;
    
    map<string,int> nodeListMap; 
    
    vector<double> jointProbs;
    
    
    void readinMatrix(string fileName, char type ); 
    void readinEdges(string edgeFileName);    
    void buildNetwork();
    double calculateProbA(int p,int c);
    double inferActivation(int p, int c);
    double getCPT1valueOfANode(int Aindex, int tumorID);
    void getCPTvalueOfOneChildOfANode(int childIndex, int Aindex,int tumorID, vector<double>& AchildCPT); 
    
    
};

#endif /* DATASET_H */

 