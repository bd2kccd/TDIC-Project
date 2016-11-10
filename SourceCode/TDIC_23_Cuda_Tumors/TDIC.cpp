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
                                vector<string> containing the two 2 global driver.  
 * @param tumorGtIndices        Input file SGA indices      
 * @param tumorGeIndices        Input file DEG indices
 *  @param tumorPosteriorMatrix  a consecutive array that carries causal score for all GT-vs-GE pairs 
 * @d_gtDataMatrix      A GPU device pointer reference, in which SGA data is stored in a tumor-by-gene matrix represented by consecutive liner array 
 * @d_geDataMatrix      A GPU device pointer reference, in which DEG data is stored in a tumor-by-gene matrix represented by consecutive liner array 
 */
//void TDIC(GTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
//        string>& mapGlobDrivers, const int tumorID, const string outPath, const float v0, int *&d_gtDataMatrix, int *&d_geDataMatrix){
void TDIC(GTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        string>& mapGlobDrivers, vector<int>  tumorGtIndices, vector<int> tumorGeIndices, const float v0, vector<float>& tumorPosteriorMatrix, bool *&d_gtDataMatrix, bool *&d_geDataMatrix){
    
    // initializations 

//    string curTumorName = gtMatrix.getTumorNameById(tumorID);
   
    // Find the GTs and GEs of the current tumor
//    vector<int> tumorGtIndices, tumorGeIndices;
    
//    gtMatrix.findGeneWithOnesInTumor(tumorID, tumorGtIndices);
//    geMatrix.findGeneWithOnesInTumor(tumorID, tumorGeIndices);
    
    // Get nGT, nGE, nTumors
    int nGT = tumorGtIndices.size();
    int nGE = tumorGeIndices.size();
    int nTumors = gtMatrix.getNTumors();
    if (nTumors != geMatrix.getNTumors()) // number of tumors should agree between gtMatrix and geMatrix
    {
        cerr << "Bug: gtMatrix and geMatrix contain different number of tumors ";
    }
    
    
    //calculate lntumorPriors
    vector<float> lntumorMutPriors;
    gtMatrix.calcLnTumorPriors(tumorGtIndices, v0, lntumorMutPriors);
    
    // Find the global drivers corresponding to the
    vector<int> tumorGlobDriverIndices;
    if (!getDEGGlobDriverIndices(gtMatrix, geMatrix, mapGlobDrivers, tumorGeIndices, tumorGlobDriverIndices))
    {
        cout << "Error occurred when retrieving global drivers";
    }

     //get Mat pointer
    bool* gtDataMatrix = gtMatrix.getMatPtr();
    bool* geDataMatrix = geMatrix.getMatPtr(); 


    
//    /*Get names of mutated GTs and GEs for this tumor
//     */
//    vector<string> gtNames;
//    vector<string> geNames;
//    
//    gtMatrix.getGeneNamesByIndices(tumorGtIndices, gtNames);
//    geMatrix.getGeneNamesByIndices(tumorGeIndices, geNames);
  
//    vector<float> tumorPosteriorMatrix(nGT*nGE,0.0);
   
    
     //****cuda starts***********************************************************************
    // Prepare device variables and invoke kernel
    // d_cancerTypes, d_gtDataMatrix, d_geDataMatrix are not tumor specific. They were malloc and copy in the TDIC_Main, before tumor loop
    //***************************
    // malloc device global memory
    int *d_tumorGeIndices, *d_tumorGtIndices, *d_tumorGlobDriverIndices;
//    int *d_gtDataMatrix, *d_geDataMatrix;   //They are not tumor specified. Memory copy before Tumor loop
    float *d_lntumorMutPriors ,*d_tumorPosteriorMatrix;
    

    cudaMalloc((int**)&d_tumorGeIndices, nGE*sizeof(int));
    cudaMalloc((int**)&d_tumorGtIndices, nGT*sizeof(int));
    cudaMalloc((int**)&d_tumorGlobDriverIndices, nGE*sizeof(int));
    
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
    cudaMemcpy(d_tumorGlobDriverIndices, &tumorGlobDriverIndices.front(), nGE*sizeof(int), cudaMemcpyHostToDevice);
    
//    cudaMemcpy(d_gtDataMatrix, gtDataMatrix, numColT*numRowT*sizeof(int), cudaMemcpyHostToDevice);
//    cudaMemcpy(d_geDataMatrix, geDataMatrix, numColE*numRowE*sizeof(int), cudaMemcpyHostToDevice);    
    cudaMemcpy(d_tumorPosteriorMatrix, &tumorPosteriorMatrix.front(), nGT*nGE*sizeof(float), cudaMemcpyHostToDevice);
    
    cudaMemcpy(d_lntumorMutPriors, &lntumorMutPriors.front(), nGT*sizeof(float), cudaMemcpyHostToDevice);
  
    
    //invoke kernel at host side
    dim3 block (128);
    dim3 grid  ((nGE + block.x - 1) / block.x);
    
    TDIC_GeLoop<<<grid, block>>>(nGE, nGT, nTumors, d_tumorGeIndices, d_tumorGtIndices, d_tumorGlobDriverIndices, d_gtDataMatrix, d_geDataMatrix, d_lntumorMutPriors ,d_tumorPosteriorMatrix);
    
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
}



/**
 * This function parse the text file that list top 2 global drivers for each of 
 * DEGs observed in a DEG matrix. 
 * @param A string fileName
 * @return A boolean value indicating the success  
 */
bool parseGlobDriverDict(string fileName, map<string, string>& globDriverMap){
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
            globDriverMap.insert(std::pair<string, string>(fields.at(0), fields.at(1)));
                
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

bool getDEGGlobDriverIndices(GTMatrix& gtMat, TDIMatrix& geMat, map<string, string>& mapGlobDrivers, vector<int>& inDEGIndices, vector<int>& OutGlobDriverVec)
{
    /*
     * First we must get the names of the DEGs corresponding to the indices in "inDEGIndices".
     * Then, using these DEG names, we can access their global driver through our map "mapGlobDrivers" 
     * and push them onto 'OutGlobDriverVec'.
     */
    //cout << "Inside getDEGGlobDriver.\n";
    vector<string> inDEGNames;
    geMat.getGeneNamesByIndices(inDEGIndices, inDEGNames);

    vector<string> globalDriverNames;
    for(int i = 0; i < inDEGNames.size(); i++)
    {
        globalDriverNames.push_back(mapGlobDrivers[inDEGNames[i]]);
    }
    
    gtMat.getGeneIndicesByNames(globalDriverNames, OutGlobDriverVec);
    return true;
}



void TDIC_Load(string inputFileName, vector< vector<string> >& gtGeneNames, vector< vector<string> >& geGeneNames, vector<string>& tumorNames)
{
    std::stringstream ss;
    std::string line;
    ifstream inFileStream;   
    
    try{
        inFileStream.open(inputFileName.c_str());  
        if ( (inFileStream.rdstate() & std::ifstream::failbit ) != 0 )
        {
            std::cerr << "Error opening file when loading input files in TDICLoad load function, quit.\n";
            inFileStream.close();
            exit(EXIT_FAILURE);
        }
    }
    catch (...) { //std::ifstream::failure e
        cerr << "Fail to open file " << inputFileName;
        
    } 
    
    int rowInFile = 0;
    vector<string> tumorGtNames; 
    vector<string> tumorGeNames;
    while (getline(inFileStream, line)){
        if (!line.empty() && line[line.size() - 1] == '\r')
            line.erase(line.size() - 1);
//        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        stringstream ss (line);
        string tmp;
        rowInFile ++;
        if (rowInFile%4 == 1){ //tumor name
            getline(ss, tmp, ',');
            tumorNames.push_back(tmp);
            tumorGtNames.clear(); 
            tumorGeNames.clear();
            tumorGtNames.push_back("A0");
        }

        else if (rowInFile%4 == 2) //cancer type
            continue; //currently do nothing
//                bool has_only_digits = (tmp.find_first_not_of( "0123456789" ) == string::npos);
//                if (!has_only_digits)
//                {
//                    cerr << "Error: Input link file second line supposed to be a digit which is a cancer type/n";
//                    exit(1);
//                }
//                cancerType = atoi(tmp.c_str());
//                if (cancerType == 0)
//                    hasCanType = false;
        
        else if (rowInFile%4 == 3){ //Gt names
            while (getline(ss, tmp, ',')){
                tumorGtNames.push_back(tmp);
            }
        }

        else if (rowInFile%4 == 0){ //Ge names
            while (getline(ss, tmp, ',')){
                tumorGeNames.push_back(tmp);
            }
            
            //push tumorGtNames, tumorGeNames into gtGeneNames, geGeneNames;
            gtGeneNames.push_back(tumorGtNames);
            geGeneNames.push_back(tumorGeNames);
        }

    }
    inFileStream.close();   

}


/**
 * This function write TDIC or PanCanTDIC results to a csv file
 * @param tumorPosteriorMatrix: a vector contains the results of the TDIC which will be passed into this function to be saved to file
 * @param curTumorName: a string contains the current tumor name of the output
 * @param gtNames: a vector contains the gt gene names of the outputs
 * @param geNames: a vector contains the ge gene names of the outputs.
 * @param outPath: a string contains the file path of the output files
 */
void TDIC_Output(vector<float>& tumorPosteriorMatrix, string curTumorName, vector<string> gtNames, vector<string> geNames, string outPath){
    // save results to file
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
    int nGE = geNames.size();
    for(int i = 0; i < nGE; i++)
    {
        outFile << "," << geNames[i];
    }
    outFile << "\n";
    
    int nGT = gtNames.size();
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
    
}

/**
 * This function pass in the gene names (gt/ge) of a tumor, and get their gene indices corresponding to TDIC ge/gt matrix 
 * @param inGeneNames(in): a vector contains the gene names of a input tumor
 * @param matrixGeneNames(in): a vector contains the gene Names of TDIC gt/ge matrix
 * @param gtIndices(out): a vector contains input genes indices corresponding to TDIC gtMatrix or geMatrix
 * 
 */
void getTumorGeneIndices(vector<string>& inGeneNames, vector<string>& matrixGeneNames, vector<int>& outGeneIndices)
{   
    for(int i = 0; i < inGeneNames.size(); i++)
    {
 
        int iFound = 0;
        for(int j = 0; j < matrixGeneNames.size(); j++)
        {

            if(inGeneNames[i].compare(matrixGeneNames[j]) == 0)
            {
                outGeneIndices.push_back(j);
                iFound = 1;
                break;
            }
        }
        if (iFound == 0)
            cout << "Customer input gene "<< inGeneNames[i] << " is not found in TDIC database\n";
    }
    
    if(inGeneNames.size() != outGeneIndices.size())
    {
        cout << "Some genes in the customer input file do not exist in TDIC database. \n";
    }
}