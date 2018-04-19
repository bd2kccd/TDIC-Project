/* 
 * File:   GTMatrix.cpp
 * Author: kevin
 * 
 * Created on May 21, 2015, 9:36 AM
 */

#include <math.h>

//#include "TDIMatrix.h"
#include "GTMatrix.h"

using namespace std;

#include <algorithm>    // std::find
#include <vector>
#include <fstream> 
#include <iostream>
#include <string>
#include <string.h>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <float.h>
#include <cstring>
#include <cstdlib>

GTMatrix::GTMatrix() : TDIMatrix() {
//    allAlphaPriors = NULL;
    //this constructor won't execute.

}
GTMatrix::GTMatrix(string gtfileName, string priorfileName) 
{
    mat_prior = NULL;
    this->load(gtfileName,priorfileName);
//    calcGlobalAlphaPriors();
}

GTMatrix::~GTMatrix(){
//    if(allAlphaPriors) delete allAlphaPriors;
    if (mat_prior) delete [] mat_prior; 
}

/**
 * Create array representing prior probability for all GTs.
 */
void GTMatrix::calcGlobalAlphaPriors()
{
    
    //initialize alpha prior array
    allAlphaPriors = new float[nGenes]();
    
    /*
     * Iterate through all tumors of the GTMatrix.
     * For each tumor, find GTs that are mutated for that tumor.
     * Then calculate a normalized weight for the probability of one of these
     * GTs being the driving mutation. This weight is 1 / (number of mutated GTs).
     * Finally, update indeces in "allAlphaPriors" corresponding to these mutated
     * GTs by adding this normalized weight to the existing value in allAlphaPriors[GT].
     *
     */
    vector<int> geneIndices;
    for(int i = 0; i < nTumors; i++)
    {   
        findGeneWithOnesInTumor(i, geneIndices);
        int numMutatedGTsForTumor = geneIndices.size();
        float gtWeight =   1.0 / (float)numMutatedGTsForTumor;
        for(int j = 0; j < numMutatedGTsForTumor; j++)
        {
            allAlphaPriors[geneIndices[j]] += gtWeight;
        }
        geneIndices.clear();
    }
}
/**
 * 
 * @param gtIndx List of GT indices for which you want to calculate the log prior.
 * @param v0 
 * @param lnTumorPriorAlphas List of normalized, log prior probabilities for each of the GTs in gtIndx 
 *         (with original sequence intact).
 */
void GTMatrix::calcLnTumorPriors(vector<int>& gtIndx, const float v0, vector<float>& lnTumorPriorAlphas)
{
    float sumAlphaPrior = 0.0;
    
    //calculate sum of priors to normalize the prior findGeneWithOnesInTumor
    for(int gt = 0; gt < gtIndx.size(); gt++)
    {
        if(gtIndx[gt] == 0)
        {
            continue;
        }

        sumAlphaPrior += allAlphaPriors[gtIndx[gt]];
    }
    
    for(int gt = 0; gt < gtIndx.size(); gt++)
    {
        //push ln(v0) onto lnTumorPriorAlphas first
        if(gtIndx[gt] == 0)
        {
            lnTumorPriorAlphas.push_back(log(v0));
            continue;
        }

        //for all other GTs, push back ln(prior(GT) / sumAlphaPrior * (1 - v0))
        float result = (1 - v0) * allAlphaPriors[gtIndx[gt]] / sumAlphaPrior;
        lnTumorPriorAlphas.push_back(log(result));
    }
}

/**
 * Load GT or GE data from raw data (CSV file).
 * @param fileName Filepath for a CSV file representing a GT or GE matrix.
 */
void GTMatrix::load(string gtfileName, string priorfileName){
    
    //read in gtMatrix, and save to dynamic array mat
    std::stringstream ss;
    std::string line;
    ifstream inFileStream;   
    vector<bool*> matrixAsVec;
    int nCol = 0, nRow = 0;

    try{
        inFileStream.open(gtfileName.c_str());
        if ( (inFileStream.rdstate() & std::ifstream::failbit ) != 0 )
        {
            std::cerr << "Error opening file when loading GT matrix, quit.\n";
            inFileStream.close();
            exit(EXIT_FAILURE);
        }
    }
    catch (...) { //std::ifstream::failure e
        cerr << "Fail to open gtMatrix file " << gtfileName;
        
    } 
            
    // read in first line to figure out the number of columns and column Names;
    getline(inFileStream, line);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    
    string tmp;
    stringstream firstSS (line);

    bool firstColFlag = true;
    geneNames.push_back("A0"); 
    while (getline(firstSS, tmp, ',' )){
        if(firstColFlag)
        {
                firstColFlag = false;
                continue;
        }
        geneNames.push_back(tmp);
  
        
        //geneIndxMap.insert(std::pair<string, int>(tmp, nCol));
        nCol ++;
    }
    
    //User may accidently input GTMatrix as PanCanGTMatrix, so check the input file is actual GTMatrix
    //Check the last column to see if the column name is something like 'Cancel Type', if it is , then exit program
    // Cancel type column name can be 'Cancer Type' or 'Cancel Types' or 'Can Type' or 'Can Types' or 'Cancer Type Code' or 'CanType code'. 
    //There can be no or more space between each word. Each character in the word is not case sensitive, can be lower or upper case.
   
    string colName;
    colName = geneNames.back();
    //Convert to lower case
    transform(colName.begin(), colName.end(), colName.begin(), ::tolower);
    //remove all spaces
    colName.erase(remove(colName.begin(), colName.end(), ' '), colName.end());
    
    bool hasCanTypeCol = true;
    if (colName == "cancertype" || colName == "cancertypes" || colName == "cantype" ||colName == "cantypes" 
            || colName == "cantypecode" || colName == "cancertypecode"|| colName == "cantypescode" || colName == "cancertypescode")
    {
        //last column is cancel type not gene name
        cout << "Input PanCanGTMatrix has cancer type";
        geneNames.pop_back();
        nCol--;
//        exit(1);
    }
    else
        hasCanTypeCol = false;

    while (getline(inFileStream, line)){
        firstColFlag = true; 

        stringstream anotherss (line);
        string tmp;
        int curCol = 0;
        while (getline(anotherss, tmp, ',')){
                if(firstColFlag)
                {
                    firstColFlag = false;
                    tumorNames.push_back(tmp);
                    matrixAsVec.push_back(new bool[nCol]());
                    continue;
                }
                if (hasCanTypeCol)
                {
                    if(curCol != nCol)
                    {
                        matrixAsVec[nRow][curCol] = atoi(tmp.c_str());
                        curCol++;
                    }
                }
                else
                {
                    matrixAsVec[nRow][curCol] = atoi(tmp.c_str());
                    curCol++;
                }
        }
        nRow++;     
    }
    inFileStream.close();            
    
    // transform the vector of inter arrays to a consecutive array so that 
    // computation can be done efficiently
    int indx = 0;
    mat = new bool[(nCol+1) * nRow ]();
    for (int i = 0; i < nRow; i++)
    {
        mat[indx++] = 1;
    }
    
    for (int g = 0; g < nCol; g++){
        for (int t = 0; t < nRow; t++){
            mat[indx++] = matrixAsVec[t][g];
        }
    }
    
    nTumors = tumorNames.size();
    nGenes = geneNames.size();
    
    // free temporary matrix
    for (int i = 0; i < matrixAsVec.size(); i++) 
        delete[] matrixAsVec[i];



    /////////////////////////////////////////////////////////////    
    //read in prior Matrix, and save to dynamic array matPrior
    //matrixAsVec_p has the first column as A0 prior

    std::string line_p;
    ifstream inFileStream_p;   
    vector<float*> matrixAsVec_p;
    int nCol_p = 0, nRow_p = 0;

    try{
        inFileStream_p.open(priorfileName.c_str());
        if ( (inFileStream_p.rdstate() & std::ifstream::failbit ) != 0 )
        {
            std::cerr << "Error opening file when loading prior matrix, quit.\n";
            inFileStream_p.close();
            exit(EXIT_FAILURE);
        }
    }
    catch (...) { //std::ifstream::failure e
        cerr << "Fail to open prior file " << priorfileName;
        
    } 
            
    // read in first line of row names and do nothing
    getline(inFileStream_p, line_p);
    //It is supposed that prior matrix should have the same col number with the gt matrix
    nCol_p = nCol+1;//it includes A0
    // read in the following lines
    while (getline(inFileStream_p, line_p)){
        firstColFlag = true; 

        stringstream ss (line_p);
        string tmp;
        int curCol = 0;
        double f_prior;
        while (getline(ss, tmp, ',')){
                if(firstColFlag)
                {
                    firstColFlag = false;
                    matrixAsVec_p.push_back(new float[nCol_p]());
                    continue;
                }
                stringstream s_tmp;
                s_tmp<<tmp;
                s_tmp>>f_prior;
                matrixAsVec_p[nRow_p][curCol] = f_prior;
                curCol++;

        }
        nRow_p++;     
    }
    inFileStream_p.close();    
    
    if (nRow_p != nRow)
    {
        cerr << "piror matrix has differenct column number than Gt matrix which is not correct!";
        exit(1);
    }
    // transform the vector of inter arrays to a consecutive array so that 
    // computation can be done efficiently
    int indx_p = 0;
    mat_prior = new float[(nCol_p) * nRow ]();
    
    for (int g = 0; g < nCol_p; g++){
        for (int t = 0; t < nRow_p; t++){
                mat_prior[indx_p++] = matrixAsVec_p[t][g];
                
        }
    }
    
    // free temporary matrix
    for (int i = 0; i < matrixAsVec_p.size(); i++) 
        delete[] matrixAsVec_p[i];


}

