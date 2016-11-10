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
    void TDIC(GTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        string>& mapGlobDrivers, vector<int>  tumorGtIndices, vector<int> tumorGeIndices, const float v0, vector<float>& tumorPosteriorMatrix, 
            bool *&d_gtDataMatrix, bool *&d_geDataMatrix);
    void TDIC_Load(string inputFileName, vector< vector<string> >& gtGeneNames, vector< vector<string> >& geGeneNames, vector<string>& tumorNames);
    void getTumorGeneIndices(vector<string>& inGeneNames, vector<string>& matrixGeneNames, vector<int>& outGeneIndices);
    //bool parseGlobDriverDict(string fileName, map<string, string> globDriverMap);
    bool parseGlobDriverDict(string fileName, map<string, string>& globDriverMap);
    bool getDEGGlobDriverIndices(GTMatrix& gtMat, TDIMatrix& geMat, map<string, string>& mapGlobDrivers, vector<int>& inDEGIndices, vector<int>& OutGlobDriverVec);
    //vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
    vector<std::string> split(const std::string &s, char delim);
    void TDIC_Output(vector<float>& tumorPosteriorMatrix, string curTumorName, vector<string> gtNames, vector<string> geNames, string outPath);
    

#endif	/* TDIC_H */

