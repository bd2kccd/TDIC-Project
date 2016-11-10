/* 
 * File:   PanCanTDIC.cpp
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
 * GT-vs-GE pairs observed in a given tumor, populate a GT-by-GE score matrix .  
 * @param gtMatrix     An TDIMatrix reference, in which SGA data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param geMatrix        An TDIMatrix reference, in which DEG data of a collection of tumors are stored
 *                              in a tumor-by-gene matrix
 * @param mapGlobDrivers        A reference to a dictionary of global drivers of DEGs in the "geMatrix",
                                Each entry is keyed by the DEG gene name, and the value is a reference of a 
                                vector<string> containing the two 2 global driver.  
 * @param tumorGtIndices        Input file SGA indices      
 * @param tumorGeIndices        Input file DEG indices
 * @param v0                    a constant float
 * @param tumorPosteriorMatrix  a consecutive array that carries causal score for all GT-vs-GE pairs 
 */

void PanCanTDIC(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        string>& mapGlobDrivers, vector<int>  tumorGtIndices, vector<int> tumorGeIndices, const float v0, vector<float>& tumorPosteriorMatrix){
    
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
    
    //get nGT, nGE, nTumores
    unsigned int nGT = tumorGtIndices.size();
    unsigned int nGE = tumorGeIndices.size();
    int nTumors = gtMatrix.getNTumors();
    if (nTumors != geMatrix.getNTumors()) // number of tumors should agree between gtMatrix and geMatrix
    {
        cerr << "Bug: gtMatrix and geMatrix contain different number of tumors ";
    }
    
    int numCanType = gtMatrix.getNumCanType();
    
    vector<float> tumorPosteriorMatrixC(nGT * nGE * numCanType, 0.0);//save each cancer type tumorPosteriorMatrix

    // loop through each GE
    #pragma omp parallel for
    for(unsigned int ge = 0; ge < nGE; ge++)
    {
        //C[] saves the count numbers of the each  cancer type
        float* C = new float[numCanType]();
        unsigned int curGeIndx = tumorGeIndices[ge];
        unsigned int rowStartForGE = curGeIndx * nTumors; 
        
        // find the globDriver for this given ge   
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

            //Variables to save the count considering the cancel types
            float* CT = new float[2*numCanType]();
            float* CTE = new float[4*numCanType]();
            float* CTD = new float[4*numCanType]();
            float* CTDE = new float[8*numCanType]();
            
            int curGTIndx = tumorGtIndices[gt];
            int gtRowStart = curGTIndx * nTumors;
            
            //count CTDE
            for(int t = 0; t < nTumors; t++)
            {
                
                int tumorCanType = gtMatrix.getCanTypeByTumorId(t);
              
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
                
                //currentTumorCanType is from 1 to numCanType, the array starts from 0, so needs to minus 1
                CTDE[(tumorCanType-1)*8+tVal*4+dVal*2+eVal]++;
            }
            
            //Calculate CTE,CTD from CTDE
            for(int a=0; a<numCanType; a++)
                for(int b=0; b<=1; b++)
                    for(int c=0; c<=1; c++)
                    {
                        CTE[ a*4 + b*2 + c ] = CTDE[ a*8 + b*4 + 0*2 + c ] + CTDE[ a*8 + b*4 + 1*2 + c ];
                        CTD[ a*4 + b*2 + c ] = CTDE[ a*8 + b*4 + c*2 + 0 ] + CTDE[ a*8 + b*4 + c*2 + 1 ];
                    }
            //Calculate TDE from CTDE
            for(int a=0; a<=1; a++)
                for(int b=0; b<=1; b++)
                    for(int c=0; c<=1; c++)
                          for (int i=0; i<numCanType; i++) //there are numCanType cancer types
                             TDE[ a*4 + b*2 + c ] += CTDE[ i*8 + a*4 + b*2 + c ] ;
    
            //Calculate TD, TE from TDE
            for(int a=0; a<=1; a++)
                for(int b=0; b<=1; b++)
                {
                    TD[ a*2 + b ] = TDE[ a*4 + b*2 + 0 ] + TDE[ a*4 + b*2 + 1 ];
                    TE[ a*2 + b ] = TDE[ a*4 + 0*2 + b ] + TDE[ a*4 + 1*2 + b ];                            
                }
            
            //Calculate CT from CTE
            for(int a=0; a<numCanType; a++)
                for(int b=0; b<=1; b++)
                    CT[ a*2 + b ] = CTE[ a*4 + b*2 + 0 ] + CTE[ a*4 + b*2 + 1 ];
                            
                            
            //Calculate T from TE
            for(int a=0; a<=1 ;a++)   
                T[a] = TE[ a*2 + 0 ] + TE[ a*2 + 1 ];
                            
            //Calculate C[] from CT
            for (int i=0; i<numCanType; i++)
                C[i] = CT[ i*2 + 0 ] + CT[ i*2 + 1 ];
             
            //There is no count for T0
            T[0] = 0.0; 
            //There is no count for T0ge0, T0ge1
            TE[0]=TE[1] = 0.0;
            //There is no count for C0T0, C1T0, ....C15T0
            for (int i=0; i<numCanType; i++)
                CT[ i*2 + 0 ] =  0.0;
            //There is no count for C0T0E0, C0T0E1, C1TOE0, C1T0E1,.....C15T0E0, C15T0E1
            for (int i=0; i<numCanType; i++)
                CTE[ i*4 + 0*2 + 0 ] = CTE[ i*4 + 0*2 + 1 ] = 0.0;
              //CTE[0] = CTE[1] = CTE[4] = CTE[5] = 0.0;
              

            
           /**************************************************
            * Calculate lnData NOT considering cancer type 
            */
            float lnData = 0.0;
            
            float TFscore ;
            float DFscore ;
            float pGT1GE1 ;
            float pGT0GE1 ;
     
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
                lnData = -FLT_MAX;
            }
            
            //save lnData 
            tumorPosteriorMatrix[gt* nGE + ge] = lnData; 
            
          
            /*****************************************************************
             *Calculate lnData for each cancel type , save the lnData into lnDataC[]
             */
 
  
            float* lnDataC = new float[numCanType]();
            float* TFscoreC = new float[numCanType]();
            float* DFscoreC = new float[numCanType]();
            float* pGT1GE1C = new float[numCanType]();
            float* pGT0GE1C = new float[numCanType]();
            
            
 
            for (int i=0; i<numCanType; i++)
            {
                if(curGTIndx == 0)
                {
                    //TFscore = calcA0Fscore(T1,   T1ge1,  T1ge0,   T0,   T0ge1,  T0ge0);//Original code
                    //corresponding to      C0T1,  C0T1E1, C0T1E0, C0T0,  C0T0E1, C0T0E0
                    //corresponding to      CT[1]] CTE[3]  CTE[2]  CT[0]] CTE[1]  CTE[0]]
                    //TFscore = calcA0Fscore(T[1],  TE[3], TE[2], T[0],  TE[1], TE[0]);//this corresponding to C=0 | GPU according to Original code
                    //TFscore = calcA0Fscore(CT[3],  CTE[7], CTE[6], CT[2],  CTE[5], CTE[4]);//this corresponding to C=1 | GPU for cancer type 1 larger
                    TFscoreC[i] = calcA0Fscore( CT[1+2*i], CTE[3+4*i], CTE[2+4*i], CT[0+2*i], CTE[1+4*i], CTE[0+4*i] );
                }

                else 
                {
                    //TFscore = calcFscore( T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0 );
                    //corresponding to      CT[1]] CTE[3]  CTE[2]  CT[0]] CTE[1]  CTE[0]]
                    //TFscore = calcFscore( CT[3],  CTE[7], CTE[6], CT[2],  CTE[5], CTE[4] );//this corresponding to C=1 | GPU for cancer type 1 larger
                    TFscoreC[i] = calcFscore(  CT[1+2*i], CTE[3+4*i], CTE[2+4*i], CT[0+2*i], CTE[1+4*i], CTE[0+4*i] );
                }

                //float DFscore = calcFscore( D1,     D1ge1,    D1ge0,     D0,     D0ge1,     D0ge0 );//Original code
                //corresponding to           C1T0D1  C1T0D1E1  C1T0D1E0   C1T0D0   C1T0D0E1   C1T0D0E0
                //corresponding to           CTD[5]  CTDE[11]  CTDE[10]   CTD[4]   CTDE[9]    CTDE[8]
                //DFscore = calcFscore( TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0] );//GPU according to Original code
                //DFscore = calcFscore( CTD[5], CTDE[11], CTDE[10], CTD[4], CTDE[9], CTDE[8] );//GPU for cancer type 1 larger
                DFscoreC[i] = calcFscore( CTD[1+4*i], CTDE[3+8*i], CTDE[2+8*i], CTD[0+4*i], CTDE[1+8*i], CTDE[0+8*i] );

                lnDataC[i]= TFscoreC[i] + DFscoreC[i] + lntumorMutPriors[gt];

                //tumorPosteriorMatrix[gt * nGE + ge] = lnData;    Comment out for cancel type calculation; Not ready to save yet. Save lnData at the last step

                if(gt == 0)
                {
                    //pGT1GE1 = (ALPHANULL + T1ge1) / (ALPHANULL + ALPHANULL + T1);  //Original code
                    //pGT0GE1 = (ALPHANULL + D0ge1 + D1ge1) / (ALPHANULL + ALPHANULL + nTumors - T1); //Original code

                    //pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]); //GPU according to Original code
                    //pGT0GE1 = (ALPHANULL + TDE[1] + TDE[3]) / (ALPHANULL + ALPHANULL + nTumors - T[1]);//GPU according to Original code
                    //pGT1GE1 = (ALPHANULL + CTE[7]) / (ALPHANULL + ALPHANULL + CT[3]);//GPU for cancer type 1 larger
                    //pGT0GE1 = (ALPHANULL + CTDE[9] + CTDE[11]) / (ALPHANULL + ALPHANULL + nTumors - CT[3]);//GPU for cancer type 1 larger
                    pGT1GE1C[i] = (ALPHANULL + CTE[3+4*i]) / (ALPHANULL + ALPHANULL + CT[1+2*i]);
                    pGT0GE1C[i] = (ALPHANULL + CTDE[1+8*i] + CTDE[3+8*i]) / (ALPHANULL + ALPHANULL + nTumors - CT[1+2*i]);
                }
                else
                {
                    //pGT1GE1 = (ALPHAIJK11 + T1ge1) / (ALPHAIJK11 + ALPHAIJK10 + T1);//Original code
                    //pGT0GE1 = (ALPHAIJK01 + D0ge1 + D1ge1) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T1);   //Original code

                    //pGT1GE1 = (ALPHAIJK11 + TE[3]) / (ALPHAIJK11 + ALPHAIJK10 + T[1]);//GPU according to Original code
                    //pGT0GE1 = (ALPHAIJK01 + TDE[1] + TDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T[1]);   //GPU according to Original code     
                    //pGT1GE1 = (ALPHAIJK11 + CTE[7]) / (ALPHAIJK11 + ALPHAIJK10 + CT[3]);//GPU for cancer type 1 larger
                    //pGT0GE1 = (ALPHAIJK01 + CTDE[9] + CTDE[11]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - CT[3]); //GPU for cancer type 1 larger
                    pGT1GE1C[i] = (ALPHAIJK11 + CTE[3+4*i]) / (ALPHAIJK11 + ALPHAIJK10 + CT[1+2*i]);
                    pGT0GE1C[i] = (ALPHAIJK01 + CTDE[1+8*i] + CTDE[3+8*i]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - CT[1+2*i]);    
                }

                if(pGT1GE1C[i] <= pGT0GE1C[i])
                {
                    //tumorPosteriorMatrix[gt* nGE + ge] = -FLT_MAX;      Changed for cancel type calculation. Save lnData at the last step
                    lnDataC[i] = -FLT_MAX;
                }
                //save lnDataC[]]
                //tumorPosteriorMatrixC0[gt* nGE + ge] = lnDataC0;
                tumorPosteriorMatrixC[ (gt* nGE + ge ) + (nGT * nGE) * i ] = lnDataC[i];
                
            }
            //delete dynamic arrays
            delete[] CT;
            delete[] CTE;
            delete[] CTD;
            delete[] CTDE;
            delete[] lnDataC;
            delete[] TFscoreC;
            delete[] DFscoreC;
            delete[] pGT1GE1C;
            delete[] pGT0GE1C;
                
        }
        
        float normalizer = 0.0;
        float* normalizerC = new float[numCanType]();
        
        
        for(unsigned int gt = 0; gt < nGT; gt++)
        {
            if(gt == 0)
            {
                normalizer = tumorPosteriorMatrix[gt * nGE + ge];
                for (int i=0; i<numCanType; i++)
                    normalizerC[i] = tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i  ]; //each cancer type shifts nGT * nGE array index
             }
            else
            {
                normalizer = logSum(normalizer, tumorPosteriorMatrix[gt * nGE + ge]);
                for (int i=0; i<numCanType; i++)
                    normalizerC[i] = logSum(normalizerC[i], tumorPosteriorMatrixC[ (gt * nGE + ge) + (nGT * nGE) * i ]);
             }
        }
        
        // finished populating a column of GTs with respect to a given GE, normalize so that the column sum to 1
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            tumorPosteriorMatrix[gt * nGE + ge] = exp(tumorPosteriorMatrix[gt * nGE + ge] - normalizer);     
            for (int i=0; i<numCanType; i++)
               tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i] = exp(tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i] - normalizerC[i]);  
        }
        
 
        float normT = 0.0; //without considering cancer type
        float* normC = new float[numCanType]();
        float normCT = 0.0; //combination of  cancer type normC[]
       
        
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
 
            normT =  tumorPosteriorMatrix[gt * nGE + ge];  

            normCT = 0.0;
            for (int i=0; i<numCanType; i++)
            {
                normC[i] = tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i];
                normCT += normC[i] * C[i]/nTumors;
            }          
                
            tumorPosteriorMatrix[gt * nGE + ge] = (normCT + normT)/2;
                    
        }  
        //delete dynamic arrays
        delete[] C;
        delete[] normalizerC;
        delete[] normC;        

    }
  
}
