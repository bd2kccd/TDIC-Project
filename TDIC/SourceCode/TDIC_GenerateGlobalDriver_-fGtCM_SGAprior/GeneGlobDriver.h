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
//    void TDIC_Global(GTMatrix&, TDIMatrix&, map<string, string> & , const int, const string outPath, const float v0);
    void GeneGlobDriver(GTMatrix&, TDIMatrix&, const string outPath, const float v0);
    //bool parseGlobDriverDict(string fileName, map<string, string> globDriverMap);
//    bool parseGlobDriverDict(string fileName, map<string, string>& globDriverMap);
//    bool getDEGGlobDriverIndices(GTMatrix& gtMat, TDIMatrix& geMat, map<string, string>& mapGlobDrivers, vector<int>& inDEGIndices, vector<int>& OutGlobDriverVec);

    //vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
    vector<std::string> split(const std::string &s, char delim);
    
    float logSum(float lnx, float lny);
    float calcFscore(float gt1,  float gt1ge1, float gt1ge0, float gt0, float gt0ge1, float gt0ge0 );
    float calcA0Fscore(float gt1,  float gt1ge1, float gt1ge0, float gt0, float gt0ge1, float gt0ge0 );

    
    
    
    
    




#endif	/* TDIC_H */

