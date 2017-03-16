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
#include "TDIC_GeLoop_Cuda.h"

//#include "TDIMatrix.h"
#include "GTMatrix.h"

#ifndef TDIC_H
#define	TDIC_H

    
    //Function declarations
//    void TDIC(GTMatrix&, TDIMatrix&, map<string, string> & , const int, const string outPath, const float v0,  int*& , int*&);
    void TDIC(GTMatrix&, TDIMatrix&, map<string, string> & , const int, const string outPath, const float v0,  bool*& , bool*&);
    //bool parseGlobDriverDict(string fileName, map<string, string> globDriverMap);
    bool parseGlobDriverDict(string fileName, map<string, string>& globDriverMap);
    bool getDEGGlobDriverIndices(GTMatrix& gtMat, TDIMatrix& geMat, map<string, string>& mapGlobDrivers, vector<int>& inDEGIndices, vector<int>& OutGlobDriverVec);
    int getGlobDriver4GE(map<string, string>& mapGlobDrivers, int geId);
    //vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
    vector<std::string> split(const std::string &s, char delim);
    

#endif	/* TDIC_H */

