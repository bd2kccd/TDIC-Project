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
 * @d_gtDataMatrix      A GPU device pointer reference, in which SGA data is stored in a tumor-by-gene matrix represented by consecutive liner array 
 * @d_geDataMatrix      A GPU device pointer reference, in which DEG data is stored in a tumor-by-gene matrix represented by consecutive liner array 
 */
//void PanCanTDIC(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
//        string>& mapGlobDrivers, const int tumorID, const string outPath, const float v0, int nCanType, int  *&d_cancerTypes,  int *&d_gtDataMatrix, int *&d_geDataMatrix){
void PanCanTDIC(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        string>& mapGlobDrivers, const int tumorID, const string outPath, const float v0, int nCanType, int  *&d_cancerTypes,  bool *&d_gtDataMatrix, bool *&d_geDataMatrix){
    // initializations 
    string curTumorName = gtMatrix.getTumorNameById(tumorID);
   
    // Find the GTs and GEs of the current tumor
    vector<int> tumorGtIndices, tumorGeIndices;
    gtMatrix.findGeneWithOnesInTumor(tumorID, tumorGtIndices);
    geMatrix.findGeneWithOnesInTumor(tumorID, tumorGeIndices);
    
//    int* gtDataMatrix = gtMatrix.getMatPtr();
//    int* geDataMatrix = geMatrix.getMatPtr();
    bool* gtDataMatrix = gtMatrix.getMatPtr();
    bool* geDataMatrix = geMatrix.getMatPtr();
    
    // Find the global drivers corresponding to the
    vector<int> tumorGlobDriverIndices;
    //map <string, string> globalDriversMap;
    if (!getDEGGlobDriverIndices(gtMatrix, geMatrix, mapGlobDrivers, tumorGeIndices, tumorGlobDriverIndices))
    {
        cout << "Error occurred when retrieving global drivers";
    }

    // Allocate memory for nGT x nGE matrix
    unsigned int nGT = tumorGtIndices.size();
    unsigned int nGE = tumorGeIndices.size();
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
    

    vector<float> tumorPosteriorMatrix(nGT * nGE,0.0);

    
    //****cuda starts***********************************************************************
    // Prepare device variables and invoke kernel

    //define device variables
    int *d_tumorGeIndices, *d_tumorGtIndices, *d_tumorGlobDriverIndices;
    float *d_lntumorMutPriors ,*d_tumorPosteriorMatrix;
    
    // malloc device global memory
    cudaMalloc((int**)&d_tumorGeIndices, nGE*sizeof(int));
    cudaMalloc((int**)&d_tumorGtIndices, nGT*sizeof(int));
    cudaMalloc((int**)&d_tumorGlobDriverIndices, nGE*sizeof(int));
         
    cudaMalloc((float**)&d_lntumorMutPriors, nGT*sizeof(float));
    cudaMalloc((float**)&d_tumorPosteriorMatrix, nGT*nGE*sizeof(float));

    
    //transfer data from host to device
    cudaMemcpy(d_tumorGeIndices, &tumorGeIndices[0], nGE*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tumorGtIndices, &tumorGtIndices[0], nGT*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_tumorGlobDriverIndices, &tumorGlobDriverIndices[0], nGE*sizeof(int), cudaMemcpyHostToDevice);
    
    cudaMemcpy(d_lntumorMutPriors, &lntumorMutPriors[0], nGT*sizeof(float), cudaMemcpyHostToDevice);

//    size_t size;
//    cudaDeviceGetLimit(&size, cudaLimitMallocHeapSize);
//    printf("Heap size found to be %d\n",(int)size);     
  
    //invoke kernel at host side
    dim3 block (128);
    dim3 grid  ((nGE + block.x - 1) / block.x);
    PanCanTDIC_GeLoop_Cuda<<<grid, block>>>(nGE, nGT, nTumors, nCanType,d_tumorGeIndices, d_tumorGtIndices, d_tumorGlobDriverIndices, d_cancerTypes, d_gtDataMatrix, d_geDataMatrix, d_lntumorMutPriors ,d_tumorPosteriorMatrix);
   
    //copy kernel result back to host side
    cudaMemcpy( &tumorPosteriorMatrix[0],d_tumorPosteriorMatrix, nGT*nGE*sizeof(float), cudaMemcpyDeviceToHost);

//    cudaDeviceSynchronize();
//    cout <<  cudaGetErrorString (cudaGetLastError() )<< "\n";
    
    // free device global memory

    cudaFree(d_tumorGeIndices);
    cudaFree(d_tumorGtIndices);
    cudaFree(d_tumorGlobDriverIndices);
    cudaFree(d_lntumorMutPriors);
    cudaFree(d_tumorPosteriorMatrix);

        



	
	
	
	
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
//    delete [] tumorPosteriorMatrixC;

}
