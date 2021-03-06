/* 
 * File:   PanCanTDIC.cpp
 *
 * 
 * 
 */

//#include "TDIC.h"


//#include "TDIMatrix.h"
//#include "GTMatrix.h"
//#include "PanCanGTMatrix.h"
#include <fstream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <float.h>
#include "PanCanTDIC.h"
//#include "TDIMatrix.h"
//#include "PanCanGTMatrix.h"
#include "PanCanTDIC_GeLoop_Cuda.h"

using namespace std;
/**
 * This function performs tumor-specific driver identification.  It calculate the causal score for all 
 * GT-vs-GE pairs observed in a given tumor, populate a GT-by-GE score matrix and output to file.  
 * @param gtMatrix     An TDIMatrix reference, in which SGA data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param geMatrix        An TDIMatrix reference, in which DEG data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param mapGlobDrivers        A reference to a dictionary of global drivers of DEGs in the "geMatrix",
                                Each entry is keyed by the DEG gene name, and the value is a reference of a 
                                vector<string> containing the two 2 global driver.  
 * @param tumorID               The tumor to be process
 * @param d_cancerTypes     A GPU device pointer to a cancer type array for all tumors
 * 
 * @param tumorPosteriorMatrix  a consecutive array that carries causal score for all GT-vs-GE pairs 
 * @d_gtDataMatrix      A GPU device pointer reference, in which SGA data is stored in a tumor-by-gene matrix represented by consecutive liner array 
 * @d_geDataMatrix      A GPU device pointer reference, in which DEG data is stored in a tumor-by-gene matrix represented by consecutive liner array 
 */
//void PanCanTDIC(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
//        string>& mapGlobDrivers, const int tumorID, const string outPath, const float v0,  int  *&d_cancerTypes,  int *&d_gtDataMatrix, int *&d_geDataMatrix){
//    // initializations 
void PanCanTDIC(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        string>& mapGlobDrivers, vector<int>  tumorGtIndices, vector<int> tumorGeIndices, const float v0, vector<float>& tumorPosteriorMatrix,  
        int realNumOfCanTypes, int  *&d_cancerTypes,  bool *&d_gtDataMatrix, bool *&d_geDataMatrix, int tumorCanType){

    //get nGT, nGE, nTumores
    unsigned int nGT = tumorGtIndices.size();
    unsigned int nGE = tumorGeIndices.size();
    int nTumors = gtMatrix.getNTumors();
    if (nTumors != geMatrix.getNTumors()) // number of tumors should agree between gtMatrix and geMatrix
    {
        cerr << "Bug: gtMatrix and geMatrix contain different number of tumors ";
    }

    
   //calculate lntumorPriors
    vector<float> lntumorMutPriors;
    gtMatrix.calcLnTumorPriors(tumorGtIndices, v0, lntumorMutPriors);

    // Find the global drivers corresponding to the ge indices
    vector<int> tumorGlobDriverIndices;
    if (!getDEGGlobDriverIndices(gtMatrix, geMatrix, mapGlobDrivers, tumorGeIndices, tumorGlobDriverIndices))
    {
        cout << "Error occurred when retrieving global drivers";
    }
 
     //get Mat pointer
    bool* gtDataMatrix = gtMatrix.getMatPtr();
    bool* geDataMatrix = geMatrix.getMatPtr(); 
    
    //We can not use realNumOfCanTypes to define a static array. So we define a constant number, it may be bigger than real cancer type numbers
    const int nCanType = 23;
//    float* tumorPosteriorMatrix = new float[nGT * nGE]();
    float* tumorPosteriorMatrixC = new float[nGT * nGE * nCanType](); //for nCanType cancer type data
    
    //****cuda starts***********************************************************************
    // Prepare device variables and invoke kernel

    //define device variables
    int *d_tumorGeIndices, *d_tumorGtIndices, *d_tumorGlobDriverIndices;
    float *d_lntumorMutPriors ,*d_tumorPosteriorMatrix, *d_tumorPosteriorMatrixC;
    
    // malloc device global memory
    cudaMalloc((int**)&d_tumorGeIndices, nGE*sizeof(int));
    cudaMalloc((int**)&d_tumorGtIndices, nGT*sizeof(int));
    cudaMalloc((int**)&d_tumorGlobDriverIndices, nGE*sizeof(int));
         
    cudaMalloc((float**)&d_lntumorMutPriors, nGT*sizeof(float));
    cudaMalloc((float**)&d_tumorPosteriorMatrix, nGT*nGE*sizeof(float));
    cudaMalloc((float**)&d_tumorPosteriorMatrixC, nGT*nGE*nCanType*sizeof(float));
    
    //transfer data from host to device
    cudaMemcpy(d_tumorGeIndices, &tumorGeIndices.front(), nGE*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tumorGtIndices, &tumorGtIndices.front(), nGT*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tumorGlobDriverIndices, &tumorGlobDriverIndices.front(), nGE*sizeof(int), cudaMemcpyHostToDevice);
    
    cudaMemcpy(d_lntumorMutPriors, &lntumorMutPriors.front(), nGT*sizeof(float), cudaMemcpyHostToDevice);
  
    //invoke kernel at host side
    dim3 block (128);
    dim3 grid  ((nGE + block.x - 1) / block.x);
    PanCanTDIC_GeLoop<<<grid, block>>>(nGE, nGT, nTumors, d_tumorGeIndices, d_tumorGtIndices, d_tumorGlobDriverIndices, 
            d_cancerTypes, d_gtDataMatrix, d_geDataMatrix, d_lntumorMutPriors ,d_tumorPosteriorMatrix, d_tumorPosteriorMatrixC, tumorCanType);

    
    //copy kernel result back to host side
    cudaMemcpy( &tumorPosteriorMatrix.front(),d_tumorPosteriorMatrix, nGT*nGE*sizeof(float), cudaMemcpyDeviceToHost);
    
    // free device global memory

    cudaFree(d_tumorGeIndices);
    cudaFree(d_tumorGtIndices);
    cudaFree(d_tumorGlobDriverIndices);
    cudaFree(d_lntumorMutPriors);
    cudaFree(d_tumorPosteriorMatrix);
    cudaFree(d_tumorPosteriorMatrixC);
        
    //***end of cuda****************************************************************************
}
