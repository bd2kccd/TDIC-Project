/* 
 * File:   TDIC.cpp
 * Author: Kevin Lu
 * 
 * Created on April 25, 2015, 6:13 PM
 */

using namespace std;
#define MAINFILE  //avoid multiple definition in TDIC.h
#include "TDIC.h"
//#include "TDIMatrix.h"
//#include "GTMatrix.h"
#include <fstream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <float.h>
#include <algorithm>

#include "TDIC.h"
//#include "TDIMatrix.h"
//#include "GTMatrix.h"
#include "TDIC_GeLoop_Cuda.h"


/**
 * This function performs tumor-specific driver identification.  It calculate the causal score for all 
 * GT-vs-GE pairs observed in a given tumor, populate a GT-by-GE score matrix and output to file.  
 * @param gtMatrix     An TDIMatrix reference, in which SGA data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param geMatrix        An TDIMatrix reference, in which DEG data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param mapGlobDrivers        A reference to a dictionary of global drivers of DEGs in the "geMatrix",
                                Each entry is keyed by the DEG gene name, and the value is a reference of a 
                                vector<string> containing a set of global drivers(exclude A0).   
 * @param tumorID               The tumor to be process
 * @d_gtDataMatrix      A GPU device pointer reference, in which SGA data is stored in a tumor-by-gene matrix represented by consecutive liner array 
 * @d_geDataMatrix      A GPU device pointer reference, in which DEG data is stored in a tumor-by-gene matrix represented by consecutive liner array 
 */
//void TDIC(GTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
//        string>& mapGlobDrivers, const int tumorID, const string outPath, const float v0, int *&d_gtDataMatrix, int *&d_geDataMatrix){

void TDIC(GTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        vector<string> >& mapGlobDrivers, const int tumorID, const string outPath, const float v0, bool *&d_gtDataMatrix, bool *&d_geDataMatrix){


    // initializations 

    string curTumorName = gtMatrix.getTumorNameById(tumorID);
   
    // Find the GTs and GEs of the current tumor
    vector<int> tumorGtIndices, tumorGeIndices;
    
    gtMatrix.findGeneWithOnesInTumor(tumorID, tumorGtIndices);
    geMatrix.findGeneWithOnesInTumor(tumorID, tumorGeIndices);
    
    bool* gtDataMatrix = gtMatrix.getMatPtr();
    bool* geDataMatrix = geMatrix.getMatPtr();

    // Find the global drivers corresponding to the
//    vector<int> tumorGlobDriverIndices;
    vector< vector<int> > tumorGlobDriverIndices;
    
    if (!getDEGGlobDriverIndices(gtMatrix, geMatrix, mapGlobDrivers, tumorGeIndices, tumorGlobDriverIndices))
    {
        cout << "Error occurred when retrieving global drivers";
    }

    // Allocate memory for nGT x nGE matrix
    int nGT = tumorGtIndices.size();
    int nGE = tumorGeIndices.size();
    int nTumors = gtMatrix.getNTumors();
    if (nTumors != geMatrix.getNTumors()) // number of tumors should agree between gtMatrix and geMatrix
    {
        cerr << "Bug: gtMatrix and geMatrix contain different number of tumors ";
    }


    cout << "Processing tumor " << curTumorName << " with " << nGT << " GAs, and " << nGE << " GEs" << "\n";
    
    /*Get names of mutated GTs and GEs for this tumor
     */
    vector<string> gtNames;
    vector<string> geNames;
    
    gtMatrix.getGeneNamesByIndices(tumorGtIndices, gtNames);
    geMatrix.getGeneNamesByIndices(tumorGeIndices, geNames);
  
    vector<float> lntumorMutPriors;
    gtMatrix.calcLnTumorPriors(tumorGtIndices, v0, lntumorMutPriors);

//    float* tumorPosteriorMatrix = new float[nGT * nGE]();
    vector<float> tumorPosteriorMatrix(nGT*nGE,0.0);
   
    
     //****cuda starts***********************************************************************
    // Prepare device variables and invoke kernel
    // d_cancerTypes, d_gtDataMatrix, d_geDataMatrix are not tumor specific. They were malloc and copy in the TDIC_Main, before tumor loop
    //***************************
    // malloc device global memory
    
    //for d_tumorPosteriorMatrix, we need to transfer the vector< vector<int> > tumorGlobDriverIndices into a liner one dimensional array
    //find the maximal size of tumorGlobDriverIndices[i].sze()
//    int maxVecSize = 0;
    int maxVecSize = 3; //consider first 3 Gts
    int numOfVec = tumorGlobDriverIndices.size();
//    for (int i=0; i<numOfVec; i++)
//    {
//        if (tumorGlobDriverIndices[i].size() > maxVecSize)
//            maxVecSize = tumorGlobDriverIndices[i].size();
//    }
    vector<int> arr_tumorGlobDriverIndices(numOfVec*maxVecSize, 0);
    int lastIndices = 0;
    for (int i=0; i<numOfVec; i++)
    {
        for (int j=0; j<maxVecSize; j++)
        {
            if (j<tumorGlobDriverIndices[i].size())
            {
                arr_tumorGlobDriverIndices[i*maxVecSize + j] = tumorGlobDriverIndices[i][j];
                lastIndices = tumorGlobDriverIndices[i][j];
            }
            else
            {
                arr_tumorGlobDriverIndices[i*maxVecSize + j] = lastIndices;
            }
        }
    }
    
    int *d_tumorGeIndices, *d_tumorGtIndices, *d_tumorGlobDriverIndices;
//    int *d_gtDataMatrix, *d_geDataMatrix;   //They are not tumor specified. Memory copy before Tumor loop
    float *d_lntumorMutPriors ,*d_tumorPosteriorMatrix;
    

    cudaMalloc((int**)&d_tumorGeIndices, nGE*sizeof(int));
    cudaMalloc((int**)&d_tumorGtIndices, nGT*sizeof(int));
    cudaMalloc((int**)&d_tumorGlobDriverIndices, numOfVec*maxVecSize*sizeof(int));
    
//    int numColT = gtMatrix.nCol;
//    int numRowT = gtMatrix.nRow;
//    cudaMalloc( (int**)&d_gtDataMatrix,numColT*numRowT*sizeof(int) );
//
//    int numColE = geMatrix.nCol;
//    int numRowE = geMatrix.nRow;
//    cudaMalloc( (int**)&d_geDataMatrix,numColE*numRowE*sizeof(int) );
     
    cudaMalloc((float**)&d_lntumorMutPriors, nGT*sizeof(float));
    cudaMalloc((float**)&d_tumorPosteriorMatrix, nGT*nGE*sizeof(float));
    
    
    //transfer data from host to device

    cudaMemcpy(d_tumorGeIndices, &tumorGeIndices.front(), nGE*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tumorGtIndices, &tumorGtIndices.front(), nGT*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tumorGlobDriverIndices, &arr_tumorGlobDriverIndices.front(), numOfVec*maxVecSize*sizeof(int), cudaMemcpyHostToDevice);
    
//    cudaMemcpy(d_gtDataMatrix, gtDataMatrix, numColT*numRowT*sizeof(int), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_geDataMatrix, geDataMatrix, numColE*numRowE*sizeof(int), cudaMemcpyHostToDevice);    
    cudaMemcpy(d_tumorPosteriorMatrix, &tumorPosteriorMatrix.front(), nGT*nGE*sizeof(float), cudaMemcpyHostToDevice);
    
    cudaMemcpy(d_lntumorMutPriors, &lntumorMutPriors.front(), nGT*sizeof(float), cudaMemcpyHostToDevice);
  
//    delete[] arr_tumorGlobDriverIndices;
    
    //invoke kernel at host side
    dim3 block (128);
    dim3 grid  ((nGE + block.x - 1) / block.x);
    
    TDIC_GeLoop<<<grid, block>>>(nGE, nGT, nTumors, d_tumorGeIndices, d_tumorGtIndices, d_tumorGlobDriverIndices, maxVecSize, d_gtDataMatrix, d_geDataMatrix, d_lntumorMutPriors ,d_tumorPosteriorMatrix);
    
    //copy kernel result back to host side
    cudaMemcpy( &tumorPosteriorMatrix.front(),d_tumorPosteriorMatrix, nGT*nGE*sizeof(float), cudaMemcpyDeviceToHost);
    
    
    // free device global memory

    cudaFree(d_tumorGeIndices);
    cudaFree(d_tumorGtIndices);
    cudaFree(d_tumorGlobDriverIndices);
    cudaFree(d_lntumorMutPriors);
    cudaFree(d_tumorPosteriorMatrix);
   
//    cudaFree(d_gtDataMatrix);
//    cudaFree(d_geDataMatrix);
        
    //***end of cuda****************************************************************************

    // save results to file
    //vector<string> gtNames = gtMatrix.getGeneNames();

    string outFileName;
    if (*outPath.end() != '/')
    {
        outFileName = outPath + "/" +  curTumorName + ".csv";
    }
    else
    {
        outFileName = outPath + curTumorName + ".csv";
    }
        

    //ofstream file;
    ofstream outFile;
    try
    {
        outFile.open(outFileName.c_str());
    }
    catch(ofstream::failure e)
    {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
    }
    
    //start writing CSV representation of TDIMatrix    
    //write column headers
    for(int i = 0; i < nGE; i++)
    {
        outFile << "," << geNames[i];
    }
    outFile << "\n";
    
    for(int i = 0; i < nGT; i++)
    {
        outFile << gtNames[i];
        for(int j = 0; j < nGE; j++)
        {
            outFile << "," << tumorPosteriorMatrix[i * nGE + j];
        }
        outFile << "\n";
    }    
    outFile.close();

//    delete [] tumorPosteriorMatrix;
}




/**
 * This function parse the text file that list top 2 global drivers for each of 
 * DEGs observed in a DEG matrix. 
 * @param A string fileName
 * @return A boolean value indicating the success  
 */
bool parseGlobDriverDict(string fileName, map<string, vector<string> >& globDriverMap){
    ifstream inFileStream;
    string line;
    vector<string> fields;
    
   
    try {
        inFileStream.open(fileName.c_str()); 
 
        while(getline(inFileStream, line))
        {   
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end()); //added 4/14/16
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); //added 4/14/16
            fields = split(line, ',');
            string key = fields.at(0);
            for(int i=1; i<fields.size(); i++)
            {
                if (i > 6)//only consider first 3 gts
                    break;
                if (!isFloat(fields.at(i)) && fields.at(i) != "A0")
                    globDriverMap[key].push_back(fields.at(i));
            }
            if (globDriverMap[key].size() == 0)
                globDriverMap[key].push_back("A0");
//            globDriverMap.insert(std::pair<string, vector<string>>(fields.at(0), fields.at(1)));
        }
        inFileStream.close();
    }
    catch (ifstream::failure e) {
        cout << "Fail to open file " << fileName;
        return false;
    } 
    return true;
}


/**
 * Split a string by a given delimiter.
 * @param s String to be split.
 * @param delim Single character delimiter by which to split the string.
 * @param elems List containing each split substring of 's' split by 'delim'.
 * @return 
 */
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

/**
 * This split function calls '&split'. User calls this function.
 * @param s String to be split by 'delim'.
 * @param delim Character delimiter to split the string 's'.
 * @return List of substrings resulting from the split.
 */
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


/**
 * 
 * @param gtMat GT matr
 * @param geMat
 * @param mapGlobDrivers
 * @param inDEGIndices
 * @param OutGlobDriverVec
 * @return 
 */

bool getDEGGlobDriverIndices(GTMatrix& gtMat, TDIMatrix& geMat, map<string, vector<string> >& mapGlobDrivers, vector<int>& inDEGIndices, vector< vector<int> >& OutGlobDriverIndices)
{
    /*
     * First we must get the names of the DEGs corresponding to the indices in "inDEGIndices".
     * Then, using these DEG names, we can access their global driver through our map "mapGlobDrivers" 
     * and push them onto 'OutGlobDriverVec'.
     */
    //cout << "Inside getDEGGlobDriver.\n";
//    vector<string> inDEGNames;
//    geMat.getGeneNamesByIndices(inDEGIndices, inDEGNames);
//
//    vector<string> globalDriverNames;
//    for(int i = 0; i < inDEGNames.size(); i++)
//    {
//        globalDriverNames.push_back(mapGlobDrivers[inDEGNames[i]]);
//    }
//    
//    gtMat.getGeneIndicesByNames(globalDriverNames, OutGlobDriverVec);
//    return true;
    
    vector<string> inDEGNames;
    geMat.getGeneNamesByIndices(inDEGIndices, inDEGNames);

    
    for(int i = 0; i < inDEGNames.size(); i++)
    {
        vector<int> globDriverIndiecsForEachGe;
        gtMat.getGeneIndicesByNames(mapGlobDrivers[inDEGNames[i]], globDriverIndiecsForEachGe);
        OutGlobDriverIndices.push_back(globDriverIndiecsForEachGe);
    }
    
    return true;    
}

int getGlobDriver4GE(map<string, string>& mapGlobDrivers, int geId)
{
    return 0;
}


///********** logSum *********************************************************/ 
///**
// * Evaluate Ln(x + y)
// * @param lnx ln(x)
// * @param lny ln(y)
// * @return ln(x + y)
// */
//__device__ float logSum(float lnx, float lny){
//    float maxExp = -4950.0;
//
//    if(lny > lnx){                
//        float tmp = lnx;
//        lnx = lny;
//        lny = tmp;
//    }
//
//    float lnyMinusLnX = lny - lnx;
//    float lnXplusLnY;
//
//    if(lnyMinusLnX < maxExp)
//        lnXplusLnY = lnx;
//    else
//        lnXplusLnY = log(1 + exp(lnyMinusLnX)) + lnx;
//
//    return (lnXplusLnY); 
//}
//
//
///***************** calcSingleGtFscore  **************************************/
///**
// * 
// * @param gt1 
// * @param gt1ge1
// * @param gt1ge0
// * @param gt0
// * @param gt0ge1
// * @param gt0ge0
// * @return 
// */
//__device__ float calcFscore(float gt1,  float gt1ge1, float gt1ge0, 
//    float gt0, float gt0ge1, float gt0ge0 )
//{
//    // Calculation of Equation 7    
//    float glnNi0 = lgamma(ALPHAIJK00 + ALPHAIJK01) - lgamma(gt0 + ALPHAIJK00 + ALPHAIJK01);
//    float glnNi1 = lgamma(ALPHAIJK10 + ALPHAIJK11) - lgamma(gt1 + ALPHAIJK10 + ALPHAIJK11);
//
//    float fscore = glnNi0 + glnNi1;   
//    fscore += lgamma(gt0ge0 + ALPHAIJK00) - lgamma(ALPHAIJK00);
//    fscore += lgamma(gt0ge1 + ALPHAIJK01) - lgamma(ALPHAIJK01);
//    fscore += lgamma(gt1ge0 + ALPHAIJK10) - lgamma(ALPHAIJK10);
//    fscore += lgamma(gt1ge1 + ALPHAIJK11) - lgamma(ALPHAIJK11);
//
//    return (fscore);
//}
//
//
///***************** calcSingleGtFscore  **************************************/
///**
// * 
// * @param gt1
// * @param gt1ge1
// * @param gt1ge0
// * @param gt0
// * @param gt0ge1
// * @param gt0ge0
// * @return 
// */
//__device__ float calcA0Fscore(float gt1,  float gt1ge1, float gt1ge0, 
//    float gt0, float gt0ge1, float gt0ge0 )
//{
//
//    // Calculation of Equation 7    
//    float glnNi0 = lgamma( ALPHANULL + ALPHANULL) - lgamma(gt0 + ALPHANULL + ALPHANULL);
//    float glnNi1 = lgamma(ALPHAIJK10 + ALPHANULL) - lgamma(gt1 + ALPHANULL + ALPHANULL);
//
//    float fscore = glnNi0 + glnNi1;   
//    fscore += lgamma(gt0ge0 + ALPHANULL) - lgamma(ALPHANULL);
//    fscore += lgamma(gt0ge1 + ALPHANULL) - lgamma(ALPHANULL);
//    fscore += lgamma(gt1ge0 + ALPHANULL) - lgamma(ALPHANULL);
//    fscore += lgamma(gt1ge1 + ALPHANULL) - lgamma(ALPHANULL);
//
//    return (fscore);
//}


bool isFloat( string myString ) {
    std::istringstream iss(myString);
    float f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail(); 
}