/* 
 * File:   GTMatrix.cpp
 * Author: kevin
 * 
 * Created on May 21, 2015, 9:36 AM
 */

#include <math.h>

#include "TDIMatrix.h"
#include "GTMatrix.h"

using namespace std;

#include <algorithm>    // std::find
#include <vector>
#include <fstream> 
#include <iostream>
#include <string>
#include <string.h>
#include <sstream>


GTMatrix::GTMatrix() : TDIMatrix() {
    allAlphaPriors = NULL;

}
GTMatrix::GTMatrix(string fileName) 
{
    this->load(fileName);
    calcGlobalAlphaPriors();
}

GTMatrix::~GTMatrix(){
    if(allAlphaPriors) delete allAlphaPriors;
}

/**
 * Create array representing prior probability for all GTs.
 */
void GTMatrix::calcGlobalAlphaPriors()
{
    //initialize alpha prior array
    allAlphaPriors = new float[nGenes];
    
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
    
    //calculate sum of priors to normalize the prior probability
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
void GTMatrix::load(string fileName){
    std::stringstream ss;
    std::string line;
    ifstream inFileStream;   
    vector<int*> matrixAsVec;
    int nCol = 0, nRow = 0;

    try{
        inFileStream.open(fileName.c_str());
        if ( (inFileStream.rdstate() & std::ifstream::failbit ) != 0 )
        {
            std::cerr << "Error opening file when loading TDI matrix, quit.\n";
            inFileStream.close();
            exit(EXIT_FAILURE);
        }
    }
    catch (...) { //std::ifstream::failure e
        cerr << "Fail to open file " << fileName;
        
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
                    matrixAsVec.push_back(new int[nCol]);
                    continue;
                }

                matrixAsVec[nRow][curCol] = atoi(tmp.c_str());
                curCol++;
        }
        nRow++;     
    }
    inFileStream.close();            
    

    
    // transform the vector of inter arrays to a consecutive array so that 
    // computation can be done efficiently
    int indx = 0;
    mat = new int[(nCol+1) * nRow ];
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
}

