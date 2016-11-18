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
    void TDIC(GTMatrix&, TDIMatrix&, map<string, vector<string> >& , const int, const string outPath, const float v0,  bool*& , bool*&);
    //bool parseGlobDriverDict(string fileName, map<string, string> globDriverMap);
    bool parseGlobDriverDict(string fileName, map<string, vector<string> >& globDriverMap);
    bool getDEGGlobDriverIndices(GTMatrix& gtMat, TDIMatrix& geMat, map<string, vector<string> >& mapGlobDrivers, vector<int>& inDEGIndices, vector< vector<int> >& OutGlobDriverIndices);
    int getGlobDriver4GE(map<string, string>& mapGlobDrivers, int geId);
    //vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
    vector<std::string> split(const std::string &s, char delim);
    bool isFloat( string myString );

#endif	/* TDIC_H */

