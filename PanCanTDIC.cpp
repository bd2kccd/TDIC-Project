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

    float* tumorPosteriorMatrix = new float[nGT * nGE]();
    float* tumorPosteriorMatrixC0 = new float[nGT * nGE]();
    float* tumorPosteriorMatrixC1 = new float[nGT * nGE]();
    
    
    // omp_set_dynamic(0);     // Explicitly disable dynamic teams
    // omp_set_num_threads(4);
    // int nthreads = omp_get_num_threads();
    // cerr << "Number of threads available for me: " << nthreads << "\n"; 
    // loop through each GE
    for(unsigned int ge = 0; ge < nGE; ge++)
    {

        //nC1 is the number of tumors of the same cancel type
        float nC1 = 0.0;
        unsigned int curGeIndx = tumorGeIndices[ge];
        unsigned int rowStartForGE = curGeIndx * nTumors; 
        
        // find the globDriver for this give ge   
        unsigned int curGDriverIndx = tumorGlobDriverIndices[ge]; //curGDriverIndx is found by ge indx        
        unsigned int rowStartForGlobDriver = curGDriverIndx * nTumors;
        
        // loop through each GT in the tumor
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            
            //Variables to save the count NOT considering the cancel type
            float T[2] = {0.0};
            float TE[4] = {0.0};
            float TD[4] = {0.0};
            float TDE[8] = {0.0};

            //Variables to save the count considering the cancel type
            float CT[4] = {0.0};
            float CTE[8] = {0.0};
            float CTD[8] = {0.0};
            float CTDE[16] = {0.0};
            
            
            int curGTIndx = tumorGtIndices[gt];

            int gtRowStart = curGTIndx * nTumors;
            
            for(int t = 0; t < nTumors; t++)
            {
                
                int compCanType = gtMatrix.getCanTypeByTumorId(t);
                int theSameCanType = (tumorCanType == compCanType);//if cancer type is the same, the value = 1
                
                int tVal = gtDataMatrix[gtRowStart + t];
                int eVal = geDataMatrix[rowStartForGE + t];
                int dVal = gtDataMatrix[rowStartForGlobDriver + t];
                
//                //Only need to count CTDE, all other variables can be calculated from CTDE
//                
//                //count while NOT considering cancel type
//                T[tVal]++;
//                TE[tVal*2+eVal]++;
//                TD[tVal*2+dVal]++;
//                TDE[tVal*4+dVal*2+eVal]++;
//                
//                //count considering cancel type                
//                CT[theSameCanType*2+tVal]++;
//                CTE[theSameCanType*4+tVal*2+eVal]++;
//                CTD[theSameCanType*4+tVal*2+dVal]++;
                
                
                CTDE[theSameCanType*8+tVal*4+dVal*2+eVal]++;
            }
            
            //Calculate CTE,CTD,TDE from CTDE
            for(int a=0; a<=1; a++)
                for(int b=0; b<=1; b++)
                    for(int c=0; c<=1; c++)
                    {
                        CTE[ a*4 + b*2 + c ] = CTDE[ a*8 + b*4 + 0*2 + c ] + CTDE[ a*8 + b*4 + 1*2 + c ];
                        CTD[ a*4 + b*2 + c ] = CTDE[ a*8 + b*4 + c*2 + 0 ] + CTDE[ a*8 + b*4 + c*2 + 1 ];
                        TDE[ a*4 + b*2 + c ] = CTDE[ 0*8 + a*4 + b*2 + c ] + CTDE[ 1*8 + a*4 + b*2 + c ];
                    }
                
            //Calculate TD, TE from TDE
            for(int a=0; a<=1; a++)
                for(int b=0; b<=1; b++)
                {
                    TD[ a*2 + b ] = TDE[ a*4 + b*2 + 0 ] + TDE[ a*4 + b*2 + 1 ];
                    TE[ a*2 + b ] = TDE[ a*4 + 0*2 + b ] + TDE[ a*4 + 1*2 + b ];                            
                }
            
            //Calculate CT from CTE
            for(int a=0; a<=1; a++)
                for(int b=0; b<=1; b++)
                    CT[ a*2 + b ] = CTE[ a*4 + b*2 + 0 ] + CTE[ a*4 + b*2 + 1 ];
                            
                            
            //Calculate T from TE
            for(int a=0; a<=1 ;a++)   
                T[a] = TE[ a*2 + 0 ] + TE[ a*2 + 1 ];
                            
            //nC1 is the number of tumors of the same cancel type
            nC1 = CT[2]+CT[3];  // nC1= C1T0 + C1T1=CT[2]+CT[3]
             
            //There is no count for T0
            T[0] = 0.0; 
            //There is no count for T0ge0, T0ge1
            TE[0]=TE[1] = 0.0;
            //There is no count for C0T0, C1T0
            CT[0] = CT[2] = 0.0;
            //There is no count for C0T0E0, C0T0E1 C1TOE0 C1T0E1
            CTE[0] = CTE[1] = CTE[4] = CTE[5] = 0.0;
              
            float TFscore ;
            float DFscore ;
            float pGT1GE1 ;
            float pGT0GE1 ;
            
           /*
            * Calculate lnData NOT considering cancer type 
            */
            
            float lnData = 0.0;
            
     
            if(curGTIndx == 0)
            {
                //TFscore = calcA0Fscore(T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0);
                TFscore = calcA0Fscore(T[1],  TE[3], TE[2], T[0],  TE[1], TE[0]);
            }
            else 
            {
                //TFscore = calcFscore( T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0 );
                TFscore = calcFscore( T[1],  TE[3], TE[2], T[0],  TE[1], TE[0] );
            }

            //float DFscore = calcFscore( D1, D1ge1, D1ge0, D0, D0ge1, D0ge0 );
            DFscore = calcFscore( TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0] );

            lnData = TFscore + DFscore + lntumorMutPriors[gt];

            //tumorPosteriorMatrix[gt * nGE + ge] = lnData;    Comment out for cancel type calculation; Not ready to save yet. Save lnData at the last step
            if(gt == 0)
            {
                //pGT1GE1 = (ALPHANULL + T1ge1) / (ALPHANULL + ALPHANULL + T1);
                //pGT0GE1 = (ALPHANULL + D0ge1 + D1ge1) / (ALPHANULL + ALPHANULL + nTumors - T1);
                pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]);
                pGT0GE1 = (ALPHANULL + TDE[1] + TDE[3]) / (ALPHANULL + ALPHANULL + nTumors - T[1]);
            }
            else
            {
                //pGT1GE1 = (ALPHAIJK11 + T1ge1) / (ALPHAIJK11 + ALPHAIJK10 + T1);
                //pGT0GE1 = (ALPHAIJK01 + D0ge1 + D1ge1) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T1);       
                pGT1GE1 = (ALPHAIJK11 + TE[3]) / (ALPHAIJK11 + ALPHAIJK10 + T[1]);
                pGT0GE1 = (ALPHAIJK01 + TDE[1] + TDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T[1]);                      
            }

            if(pGT1GE1 <= pGT0GE1)
            {
                //tumorPosteriorMatrix[gt* nGE + ge] = -FLT_MAX;      Changed for cancel type calculation. Save lnData at the last step
                lnData = -FLT_MAX;
            }
            
           /*
            *Calculate lnData while cancel types are the same, save the lnData into lnDataC1
            */
 
            float lnDataC1 = 0.0;
            if(curGTIndx == 0)
            {
                //TFscore = calcA0Fscore(T1,   T1ge1,  T1ge0,   T0,   T0ge1,  T0ge0);
                //corresponding to      C1T1,  C1T1E1, C1T1E0, C1T0,  C1T0E1, C1T0E0
                //corresponding to      CT[3]] CTE[7]  CTE[6]  CT[2]] CTE[5]  CTE[4]]
                //TFscore = calcA0Fscore(T[1],  TE[3], TE[2], T[0],  TE[1], TE[0]);
                TFscore = calcA0Fscore(CT[3],  CTE[7], CTE[6], CT[2],  CTE[5], CTE[4]);
            }
            else 
            {
                //TFscore = calcFscore( T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0 );
                TFscore = calcFscore( CT[3],  CTE[7], CTE[6], CT[2],  CTE[5], CTE[4] );
            }

            //float DFscore = calcFscore( D1,     D1ge1,    D1ge0,     D0,     D0ge1,     D0ge0 );
            //corresponding to           C1T0D1  C1T0D1E1  C1T0D1E0   C1T0D0   C1T0D0E1   C1T0D0E0
            //corresponding to           CTD[5]  CTDE[11]  CTDE[10]   CTD[4]   CTDE[9]    CTDE[8]
            //DFscore = calcFscore( TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0] );
            DFscore = calcFscore( CTD[5], CTDE[11], CTDE[10], CTD[4], CTDE[9], CTDE[8] );
            
            lnDataC1 = TFscore + DFscore + lntumorMutPriors[gt];

            //tumorPosteriorMatrix[gt * nGE + ge] = lnData;    Comment out for cancel type calculation; Not ready to save yet. Save lnData at the last step
            
            if(gt == 0)
            {
                //pGT1GE1 = (ALPHANULL + T1ge1) / (ALPHANULL + ALPHANULL + T1);
                //pGT0GE1 = (ALPHANULL + D0ge1 + D1ge1) / (ALPHANULL + ALPHANULL + nTumors - T1);
                
                //pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]);
                //pGT0GE1 = (ALPHANULL + TDE[1] + TDE[3]) / (ALPHANULL + ALPHANULL + nTumors - T[1]);
                pGT1GE1 = (ALPHANULL + CTE[7]) / (ALPHANULL + ALPHANULL + CT[3]);
                pGT0GE1 = (ALPHANULL + CTDE[9] + CTDE[11]) / (ALPHANULL + ALPHANULL + nTumors - CT[3]);
            }
            else
            {
                //pGT1GE1 = (ALPHAIJK11 + T1ge1) / (ALPHAIJK11 + ALPHAIJK10 + T1);
                //pGT0GE1 = (ALPHAIJK01 + D0ge1 + D1ge1) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T1);   
                
                //pGT1GE1 = (ALPHAIJK11 + TE[3]) / (ALPHAIJK11 + ALPHAIJK10 + T[1]);
                //pGT0GE1 = (ALPHAIJK01 + TDE[1] + TDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T[1]);         
                pGT1GE1 = (ALPHAIJK11 + CTE[7]) / (ALPHAIJK11 + ALPHAIJK10 + CT[3]);
                pGT0GE1 = (ALPHAIJK01 + CTDE[9] + CTDE[11]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - CT[3]);         
            }
         
            if(pGT1GE1 <= pGT0GE1)
            {
                //tumorPosteriorMatrix[gt* nGE + ge] = -FLT_MAX;      Changed for cancel type calculation. Save lnData at the last step
                lnDataC1 = -FLT_MAX;
            }
            
            
            /*
            *Calculate lnData while cancel types are NOT the same, save the lnData into lnDataC0
            */
 
            float lnDataC0 = 0.0;
            if(curGTIndx == 0)
            {
                //TFscore = calcA0Fscore(T1,   T1ge1,  T1ge0,   T0,   T0ge1,  T0ge0);
                //corresponding to      C0T1,  C0T1E1, C0T1E0, C0T0,  C0T0E1, C0T0E0
                //corresponding to      CT[1]] CTE[3]  CTE[2]  CT[0]] CTE[1]  CTE[0]]
                //TFscore = calcA0Fscore(T[1],  TE[3], TE[2], T[0],  TE[1], TE[0]);
                TFscore = calcA0Fscore(CT[1],  CTE[3], CTE[2], CT[0],  CTE[1], CTE[0]);
            }
            else 
            {
                //TFscore = calcFscore( T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0 );
                TFscore = calcFscore( CT[1],  CTE[3], CTE[2], CT[0],  CTE[1], CTE[0] );
            }

            //float DFscore = calcFscore( D1,     D1ge1,    D1ge0,     D0,     D0ge1,     D0ge0 );
            //corresponding to           C0T0D1  C0T0D1E1  C0T0D1E0   C0T0D0   C0T0D0E1   C0T0D0E0
            //corresponding to           CTD[1]  CTDE[3]  CTDE[2]   CTD[0]   CTDE[1]    CTDE[0]
            //DFscore = calcFscore( TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0] );
            DFscore = calcFscore( CTD[1], CTDE[3], CTDE[2], CTD[0], CTDE[1], CTDE[0] );
            
                       
            
            lnDataC0 = TFscore + DFscore + lntumorMutPriors[gt];

              //tumorPosteriorMatrix[gt * nGE + ge] = lnData;    Comment out for cancel type calculation; Not ready to save yet. Save lnData at the last step
            if(gt == 0)
            {
                //pGT1GE1 = (ALPHANULL + T1ge1) / (ALPHANULL + ALPHANULL + T1);
                //pGT0GE1 = (ALPHANULL + D0ge1 + D1ge1) / (ALPHANULL + ALPHANULL + nTumors - T1);
            
                //pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]);
                //pGT0GE1 = (ALPHANULL + TDE[1] + TDE[3]) / (ALPHANULL + ALPHANULL + nTumors - T[1]);
                pGT1GE1 = (ALPHANULL + CTE[3]) / (ALPHANULL + ALPHANULL + CT[1]);
                pGT0GE1 = (ALPHANULL + CTDE[1] + CTDE[3]) / (ALPHANULL + ALPHANULL + nTumors - CT[1]);
            }
            else
            {
                //pGT1GE1 = (ALPHAIJK11 + T1ge1) / (ALPHAIJK11 + ALPHAIJK10 + T1);
                //pGT0GE1 = (ALPHAIJK01 + D0ge1 + D1ge1) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T1);   
              
                //pGT1GE1 = (ALPHAIJK11 + TE[3]) / (ALPHAIJK11 + ALPHAIJK10 + T[1]);
                //pGT0GE1 = (ALPHAIJK01 + TDE[1] + TDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T[1]);         
                pGT1GE1 = (ALPHAIJK11 + CTE[3]) / (ALPHAIJK11 + ALPHAIJK10 + CT[1]);
                pGT0GE1 = (ALPHAIJK01 + CTDE[1] + CTDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - CT[1]);         
            }
         
            if(pGT1GE1 <= pGT0GE1)
            {
                //tumorPosteriorMatrix[gt* nGE + ge] = -FLT_MAX;      Changed for cancel type calculation. Save lnData at the last step
                lnDataC0 = -FLT_MAX;
            }
            
             
            //save lnData lnDataC0 lnDataC1
            tumorPosteriorMatrix[gt* nGE + ge] = lnData; 
            tumorPosteriorMatrixC0[gt* nGE + ge] = lnDataC0;
            tumorPosteriorMatrixC1[gt* nGE + ge] = lnDataC1;

    
        }
        
        float normalizer = 0;
        float normalizerC0 = 0;
        float normalizerC1 = 0;
        
        for(unsigned int gt = 0; gt < nGT; gt++)
        {
            if(gt == 0)
            {
                normalizer = tumorPosteriorMatrix[gt * nGE + ge];
                normalizerC0 = tumorPosteriorMatrixC0[gt * nGE + ge];
                normalizerC1 = tumorPosteriorMatrixC1[gt * nGE + ge];
            }
            else
            {
                normalizer = logSum(normalizer, tumorPosteriorMatrix[gt * nGE + ge]);
                normalizerC0 = logSum(normalizerC0, tumorPosteriorMatrixC0[gt * nGE + ge]);
                normalizerC1 = logSum(normalizerC1, tumorPosteriorMatrixC1[gt * nGE + ge]);
            }
        }
        
        // finished populating a column of GTs with respect to a given GE, normalize so that the column sum to 1
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            tumorPosteriorMatrix[gt * nGE + ge] = exp(tumorPosteriorMatrix[gt * nGE + ge] - normalizer);     
            tumorPosteriorMatrixC0[gt * nGE + ge] = exp(tumorPosteriorMatrixC0[gt * nGE + ge] - normalizerC0);  
            tumorPosteriorMatrixC1[gt * nGE + ge] = exp(tumorPosteriorMatrixC1[gt * nGE + ge] - normalizerC1); 
        }
        
        float normC0;  //cancer types are not the same
        float normC1; //cancer type is the same
        float normT; //total tumors without considering cancer type
        
        float sumC0C1 = 0.0; //for test purpose
        float sumT = 0.0; //for test
        
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            normC0 = tumorPosteriorMatrixC0[gt * nGE + ge] ;
            normC1 = tumorPosteriorMatrixC1[gt * nGE + ge];
            normT =  tumorPosteriorMatrix[gt * nGE + ge];   
            
            tumorPosteriorMatrix[gt * nGE + ge] = (normC0 * (nTumors - nC1)/nTumors + normC1 * nC1/nTumors + normT)/2;
//            sumC0C1 += normC0 * (nTumors - nC1)/nTumors + normC1 * nC1/nTumors;//for test
//            sumT += normT; //for test
        }    
//        if (fabs(sumC0C1-1) > 0.001) //for test purpose
//            cout << "bug: in PanCanTDIC.Cpp/PanCanTDIC sumC0C1 != 1 difference is " << sumC0C1-1 << "\n"; //for test purpose
//        if (fabs(sumT - 1) > 0.001) //for test purpose
//            cout << "bug: in PanCanTDIC.Cpp/PanCanTDIC sumT != 1 difference is " << sumT-1 << "\n"; //for test purpose
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
    delete [] tumorPosteriorMatrixC0;
    delete [] tumorPosteriorMatrixC1;
}
