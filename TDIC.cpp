/* 
 * File:   TDIC.cpp
 * Author: Kevin Lu
 * 
 * Created on April 25, 2015, 6:13 PM
 */


#include "TDIC.h"
#include "TDIMatrix.h"
#include "GTMatrix.h"
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
 */
void TDIC(GTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        string>& mapGlobDrivers, const int tumorID, const string outPath, const float v0){
    // initializations 
 
    string curTumorName = gtMatrix.getTumorNameById(tumorID);
   
    // Find the GTs and GEs of the current tumor
    vector<int> tumorGtIndices, tumorGeIndices;
    gtMatrix.findGeneWithOnesInTumor(tumorID, tumorGtIndices);
    geMatrix.findGeneWithOnesInTumor(tumorID, tumorGeIndices);
    
    int* gtDataMatrix = gtMatrix.getMatPtr();
    int* geDataMatrix = geMatrix.getMatPtr();

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

    float* tumorPosteriorMatrix = new float[nGT * nGE];
    
    
    // omp_set_dynamic(0);     // Explicitly disable dynamic teams
    // omp_set_num_threads(4);
    // int nthreads = omp_get_num_threads();
    // cerr << "Number of threads available for me: " << nthreads << "\n"; 
    // loop through each GE
    for(unsigned int ge = 0; ge < nGE; ge++)
    {
        float normalizer = 0;
        unsigned int curGeIndx = tumorGeIndices[ge];
        unsigned int rowStartForGE = curGeIndx * nTumors; 
        
        // find the globDriver for this give ge   
        unsigned int curGDriverIndx = tumorGlobDriverIndices[ge]; //curGDriverIndx is found by ge indx        
        unsigned int rowStartForGlobDriver = curGDriverIndx * nTumors;
        
        // loop through each GT in the tumor
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            //cout << "GT iteration: " << gt << "\n";
            // statistics associated with current T and global driver only 
            float T1 = 0.0,   T1ge1 = 0.0, T1ge0 = 0.0, T0 = 0.0, T0ge1 = 0.0, T0ge0 = 0.0; 
            float D1 = 0.0,   D1ge1 = 0.0, D1ge0 = 0.0, D0 = 0.0, D0ge1 = 0.0, D0ge0 = 0.0;
            
            int curGTIndx = tumorGtIndices[gt];

            int gtRowStart = curGTIndx * nTumors;
            
            for(int t = 0; t < nTumors; t++)
            {
                //if GT = 1 at current tumor, update stats for GT = 1 
                if(gtDataMatrix[gtRowStart + t] == 1)
                {   
                    T1++;
                    if(geDataMatrix[rowStartForGE + t] == 1) //
                    {
                        T1ge1++;
                    }
                    else
                    {
                        T1ge0++;
                    }
                }

                //if GT != 1, we use the global driver and collect stats for globalDriver = 1
                else
                {
                    if(gtDataMatrix[rowStartForGlobDriver + t] == 1) 
                    {
                        D1++;
                        if(geDataMatrix[rowStartForGE + t] == 1) //
                        {
                            D1ge1++;
                        }
                        else
                        {
                            D1ge0++;
                        }
                    }
                    else
                    {
                        D0++;
                        if(geDataMatrix[rowStartForGE + t] == 1)
                        {
                            D0ge1++;
                        }
                        else
                        {
                            D0ge0++;
                        }
                    }
                } 
            }

            float TFscore;
            if(curGTIndx == 0)
            {
                TFscore = calcA0Fscore(T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0);
            }
            else 
            {
                TFscore = calcFscore( T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0 );
            }

            float DFscore = calcFscore( D1, D1ge1, D1ge0, D0, D0ge1, D0ge0 );

            float lnData = TFscore + DFscore + lntumorMutPriors[gt];

            tumorPosteriorMatrix[gt * nGE + ge] = lnData;

            float pGT1GE1, pGT0GE1;
            if(gt == 0)
            {
                pGT1GE1 = (ALPHANULL + T1ge1) / (ALPHANULL + ALPHANULL + T1);
                pGT0GE1 = (ALPHANULL + D0ge1 + D1ge1) / (ALPHANULL + ALPHANULL + nTumors - T1);
            }
            else
            {
                pGT1GE1 = (ALPHAIJK11 + T1ge1) / (ALPHAIJK11 + ALPHAIJK10 + T1);
                pGT0GE1 = (ALPHAIJK01 + D0ge1 + D1ge1) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T1);                
            }

            // if(ge == 4)
            // {
            //     cout << "GT1GE1: " << pGT1GE1 << " GT0GE1: " << pGT0GE1 << " T0ge1: " << T0ge1 << "\n";

            // }
            if(pGT1GE1 <= pGT0GE1)
            {
                tumorPosteriorMatrix[gt* nGE + ge] = -FLT_MAX;
            }
        }

        for(unsigned int gt = 0; gt < nGT; gt++)
        {
            if(gt == 0)
            {
                normalizer = tumorPosteriorMatrix[gt * nGE + ge];
            }
            else
            {
                normalizer = logSum(normalizer, tumorPosteriorMatrix[gt * nGE + ge]);
            }
        }
        
        // finished populating a column of GTs with respect to a given GE, normalize so that the column sum to 1
        for (unsigned int gt = 0; gt < nGT; gt++)
            tumorPosteriorMatrix[gt * nGE + ge] = exp(tumorPosteriorMatrix[gt * nGE + ge] - normalizer);     
    }
  

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

    delete [] tumorPosteriorMatrix;
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

int getGlobDriver4GE(map<string, string>& mapGlobDrivers, int geId)
{
    return 0;
}


/********** logSum *********************************************************/ 
/**
 * Evaluate Ln(x + y)
 * @param lnx ln(x)
 * @param lny ln(y)
 * @return ln(x + y)
 */
float logSum(float lnx, float lny){
    float maxExp = -4950.0;

    if(lny > lnx){                
        float tmp = lnx;
        lnx = lny;
        lny = tmp;
    }

    float lnyMinusLnX = lny - lnx;
    float lnXplusLnY;

    if(lnyMinusLnX < maxExp)
        lnXplusLnY = lnx;
    else
        lnXplusLnY = log(1 + exp(lnyMinusLnX)) + lnx;

    return (lnXplusLnY); 
}


/***************** calcSingleGtFscore  **************************************/
/**
 * 
 * @param gt1 
 * @param gt1ge1
 * @param gt1ge0
 * @param gt0
 * @param gt0ge1
 * @param gt0ge0
 * @return 
 */
float calcFscore(float gt1,  float gt1ge1, float gt1ge0, 
    float gt0, float gt0ge1, float gt0ge0 )
{
    // Calculation of Equation 7    
    float glnNi0 = lgamma(ALPHAIJK00 + ALPHAIJK01) - lgamma(gt0 + ALPHAIJK00 + ALPHAIJK01);
    float glnNi1 = lgamma(ALPHAIJK10 + ALPHAIJK11) - lgamma(gt1 + ALPHAIJK10 + ALPHAIJK11);

    float fscore = glnNi0 + glnNi1;   
    fscore += lgamma(gt0ge0 + ALPHAIJK00) - lgamma(ALPHAIJK00);
    fscore += lgamma(gt0ge1 + ALPHAIJK01) - lgamma(ALPHAIJK01);
    fscore += lgamma(gt1ge0 + ALPHAIJK10) - lgamma(ALPHAIJK10);
    fscore += lgamma(gt1ge1 + ALPHAIJK11) - lgamma(ALPHAIJK11);

    return (fscore);
}


/***************** calcSingleGtFscore  **************************************/
/**
 * 
 * @param gt1
 * @param gt1ge1
 * @param gt1ge0
 * @param gt0
 * @param gt0ge1
 * @param gt0ge0
 * @return 
 */
float calcA0Fscore(float gt1,  float gt1ge1, float gt1ge0, 
    float gt0, float gt0ge1, float gt0ge0 )
{

    // Calculation of Equation 7    
    float glnNi0 = lgamma( ALPHANULL + ALPHANULL) - lgamma(gt0 + ALPHANULL + ALPHANULL);
    float glnNi1 = lgamma(ALPHAIJK10 + ALPHANULL) - lgamma(gt1 + ALPHANULL + ALPHANULL);

    float fscore = glnNi0 + glnNi1;   
    fscore += lgamma(gt0ge0 + ALPHANULL) - lgamma(ALPHANULL);
    fscore += lgamma(gt0ge1 + ALPHANULL) - lgamma(ALPHANULL);
    fscore += lgamma(gt1ge0 + ALPHANULL) - lgamma(ALPHANULL);
    fscore += lgamma(gt1ge1 + ALPHANULL) - lgamma(ALPHANULL);

    return (fscore);
}

/**
 * 
 * @param gtm 
 * @param gem
 * @param nTumors
 * @param gtIndx
 * @param nGT
 * @param geIndx
 * @return 
 */
float calcFscoreMultGT(int* gtm, int* gem, const int nTumors,  const int* gtIndx, const unsigned int nGT, const int geIndx)
{
    int parentCombinationCounts[2 ^ nGT];
    int parentAndChildCombinationCounts[2 ^ (nGT + 1)];
    
    unsigned int geRowStart = geIndx * nTumors; 
    unsigned int gtRowStarts[nGT];
    for (int i = 0; i < nGT; i++)
    {
        gtRowStarts[i] = nTumors * gtIndx[i];
    }
    
    for(int i = 0; i < nTumors; i++)
    {
        int slot = 0;
        for(int j = 0; j < nGT; j++)
        {
            slot = slot * 2 + gtm[gtRowStarts[j] + i];
        }
        parentCombinationCounts[slot]++;
        slot = slot * 2 + gem[geRowStart + i];
        parentAndChildCombinationCounts[slot]++;
    }
    
    // calculate Equation 7'.3
    float res = 0;
    for (int j = 0; j < 2 ^ nGT; j++)
    {
        res += lgamma(ALPHAIJK00 + ALPHAIJK01) - lgamma(parentCombinationCounts[j] + ALPHAIJK00 + ALPHAIJK01);
    }
    
    res += lgamma(parentAndChildCombinationCounts[0] + ALPHAIJK00) - lgamma(ALPHAIJK00);
    res += lgamma(parentAndChildCombinationCounts[1] + ALPHAIJK01) - lgamma(ALPHAIJK01);
    float myAlpha[2] = {1.0, 2.0};
    for (int j = 2; j < 2 ^ (nGT + 1); j++)
    {
        res += lgamma(parentAndChildCombinationCounts[j] + myAlpha[j%2]) - lgamma(myAlpha[j%2]);
    }
    
    return res;
}

