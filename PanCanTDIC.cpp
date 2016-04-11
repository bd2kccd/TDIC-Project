/* 
 * File:   PanCanTDIC.cpp
 *
 * 
 * 
 */

//#include "TDIC.h"
#include "PanCanTDIC.h"
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
void PanCanTDIC(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        string>& mapGlobDrivers, const int tumorID, const string outPath, const float v0){
    // initializations 
 
    int tumorCanType = gtMatrix.getCanTypeByTumorId(tumorID);
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
                int compCanType = gtMatrix.getCanTypeByTumorId(t);
                int NotTheSameCanType = (tumorCanType != compCanType);
                
                
                // The following code copied from original TDIC.cpp. This needs to add cancel type options. Not complete yet
                
                
                
                
                
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
