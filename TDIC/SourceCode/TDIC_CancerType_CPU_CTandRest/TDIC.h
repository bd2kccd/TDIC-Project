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
    void TDIC(GTMatrix&, TDIMatrix&, map<string, string> & , const int, const string outPath, const double v0);
    //bool parseGlobDriverDict(string fileName, map<string, string> globDriverMap);
    bool parseGlobDriverDict(string fileName, map<string, string>& globDriverMap);
    bool getDEGGlobDriverIndices(GTMatrix& gtMat, TDIMatrix& geMat, map<string, string>& mapGlobDrivers, vector<int>& inDEGIndices, vector<int>& OutGlobDriverVec);

    //vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
    vector<std::string> split(const std::string &s, char delim);
    
    double logSum(double lnx, double lny);
    double calcFscore(double gt1,  double gt1ge1, double gt1ge0, double gt0, double gt0ge1, double gt0ge0 );
    double calcA0Fscore(double gt1,  double gt1ge1, double gt1ge0, double gt0, double gt0ge1, double gt0ge0 );

    
    
    
    
    




#endif	/* TDIC_H */

