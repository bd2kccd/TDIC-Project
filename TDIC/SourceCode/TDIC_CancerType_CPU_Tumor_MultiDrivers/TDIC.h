/* 
 * File:   TDIC.h
 * Author: Kevin Lu
 *
 * Created on April 25, 2015, 10:21 PM
 */
    
    
using namespace std;
#include <map>
#include <vector>
#include <string>

//#include "TDIMatrix.h"
#include "GTMatrix.h"

#ifndef TDIC_H
#define	TDIC_H

//#ifdef	__cplusplus
//extern "C" {
//#endif

    // define constant hyper parameters
    #define ALPHANULL 1.0

    #define ALPHAIJK00 2.0
    #define ALPHAIJK01 1.0
    #define ALPHAIJK10 1.0
    #define ALPHAIJK11 2.0


//#ifdef	__cplusplus
//}
//#endif
    
    //Function declarations
    void TDIC(GTMatrix&, TDIMatrix&, map<string, vector<string> >& , vector<int>  tumorGtIndices, vector<int> tumorGeIndices, const float v0, vector<float>& tumorPosteriorMatrix);
    void TDIC_Load(string inputFileName,  vector<string>& gtGeneNames, vector<string>& geGeneNames, string& curTumorName);
    void getTumorGeneIndices(vector<string>& inGeneNames, vector<string>& matrixGeneNames, vector<int>& outGeneIndices);
    void TDIC_Output(vector<float>& tumorPosteriorMatrix, string curTumorName, vector<string> gtNames, vector<string> geNames, string outPath);
    //bool parseGlobDriverDict(string fileName, map<string, string> globDriverMap);
    bool parseGlobDriverDict(string fileName, map<string, vector<string> >& globDriverMap);
    bool getDEGGlobDriverIndices(GTMatrix& gtMat, TDIMatrix& geMat, map<string, vector<string> >& mapGlobDrivers, vector<int>& inDEGIndices, vector< vector<int> >& OutGlobDriverIndices);

    //vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
    vector<std::string> split(const std::string &s, char delim);
    
    float logSum(float lnx, float lny);
    float calcFscore(float gt1,  float gt1ge1, float gt1ge0, float gt0, float gt0ge1, float gt0ge0 );
    float calcA0Fscore(float gt1,  float gt1ge1, float gt1ge0, float gt0, float gt0ge1, float gt0ge0 );
    bool isFloat( string myString );
    


#endif	/* TDIC_H */

