/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

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
#include "PanCanTDIC_GeLoop_Cuda.h"
//#include "TDIC_GeLoop_Cuda.h"

using namespace std;

 // gpu device constant variables
__constant__ float  ALPHANULL = 1.0;
__constant__ float  ALPHAIJK00 = 2.0;
__constant__ float  ALPHAIJK01 = 1.0;
__constant__ float  ALPHAIJK10 = 1.0;
__constant__ float  ALPHAIJK11 = 2.0;
   
   
 
/********** TDIC_GeLoop *********************************************************/ 
/**
 * TDIC at ge loop send to GPU for parallel computation
 * @param int nGE   Length of tumorGeIndieces. Note tumorGeIndices is a pointer now, not a vector, no way to use .size() to get nGE
 * @param int nGT   Length of tumorGtIndieces.
 * @param int nTumors   Number of total tumors
 * @param int *tumorGeIndices   Indices of Ge of given tumor
 * @param int *tumorGtIndices   Indices of Gt of given tumor
 * @param int *tumorGlobDriverIndices   Index of global driver gt for each ge of a given tumor
 * @param int* canTypes     Cancer type of each tumor
 * @param bool *gtDataMatrix    gt-by-tumor matrix represented as a 1-D array
 * @param bool *geDataMatrix    ge-by-tumor matrix represented as a 1-D array
 * @param float *lntumorMutPriors    Prior probability that an gt is a driver in a tumor
 * @param float *tumorPosteriorMatrix      nTumorGts-by-nTUmorGes matrix represented as a 1-D array. Store the posterior probability that a Gt causes a Ge. This param will be transfered out to CPU for outputting to files
 * @param float* tumorPosteriorMatrixC     Similiar to tumorPosterior, the size is number of cancer types * sizeof(tumorPosterior), contains the probability for every  cancer type. This probability will be combined into tumorPosterior. No data transfer in or out
 */     
__global__ void PanCanTDIC_GeLoop(int nGE, int nGT, int nTumors, int* tumorGeIndices,  int* tumorGtIndices, 
        int* tumorGlobDriverIndices, int* canTypes, bool* gtDataMatrix, bool* geDataMatrix, float* lntumorMutPriors, 
        float* tumorPosteriorMatrix, float* tumorPosteriorMatrixC, int tumorCanType)
{
    int ge = threadIdx.x +  blockDim.x * blockIdx.x;//pass threadId to ge variable
    if (ge<nGE)
    {

        //C[] saves the count numbers of the each  cancer type
        float C[2] = {0.0};
        int curGeIndx = tumorGeIndices[ge];
        int rowStartForGE = curGeIndx * nTumors; 
        
        // find the globDriver for this given ge   
        int curGDriverIndx = tumorGlobDriverIndices[ge]; //curGDriverIndx is found by ge indx        
        int rowStartForGlobDriver = curGDriverIndx * nTumors;
        
        // loop through each GT in the tumor
        for (int gt = 0; gt < nGT; gt++)
        {
            
//            //Variables to save the count NOT considering the cancel type
//            float T[2] = {0.0};
//            float TE[4] = {0.0};
//            float TD[4] = {0.0};
//            float TDE[8] = {0.0};

            //Variables to save the count considering  cancel types
            float CT[4] = {0.0};
            float CTE[8] = {0.0};
            float CTD[8] = {0.0};
            float CTDE[16] = {0.0}; 
            
            
            int curGTIndx = tumorGtIndices[gt];

            int gtRowStart = curGTIndx * nTumors;
            
            //count CTDE
            for(int t = 0; t < nTumors; t++)
            {
                int compCanType = canTypes[t];
                int theSameCanType = (tumorCanType == compCanType);
                
//                int tumorCanType = gtMatrix.getCanTypeByTumorId(t);
//                if (tumorCanType < 1)
//                {
//                   cerr << "\nbug: In PanCanTDIC.Cpp:  Cancer type of tumor "<< geMatrix.getTumorNameById(t)<< " is " << tumorCanType <<" Which is not valid\n";
//                   exit(1);
//                }                          
                
                
                int tVal = gtDataMatrix[gtRowStart + t];
                int eVal = geDataMatrix[rowStartForGE + t];
                int dVal = gtDataMatrix[rowStartForGlobDriver + t];
                
                //Only need to count CTDE, all other variables can be calculated from CTDE
                
                //theSameCanType is 1 for count of cancer type the same with the tumor
                CTDE[theSameCanType*8+tVal*4+dVal*2+eVal]++;
            }
            
            //Calculate CTE,CTD from CTDE
            for(int a=0; a<=1; a++)
                for(int b=0; b<=1; b++)
                    for(int c=0; c<=1; c++)
                    {
                        CTE[ a*4 + b*2 + c ] = CTDE[ a*8 + b*4 + 0*2 + c ] + CTDE[ a*8 + b*4 + 1*2 + c ];
                        CTD[ a*4 + b*2 + c ] = CTDE[ a*8 + b*4 + c*2 + 0 ] + CTDE[ a*8 + b*4 + c*2 + 1 ];
                    }
//            //Calculate TDE from CTDE
//            for(int a=0; a<=1; a++)
//                for(int b=0; b<=1; b++)
//                    for(int c=0; c<=1; c++)
//                          for (int i=0; i<=1; i++) 
//                             TDE[ a*4 + b*2 + c ] += CTDE[ i*8 + a*4 + b*2 + c ] ;
    
//            //Calculate TD, TE from TDE
//            for(int a=0; a<=1; a++)
//                for(int b=0; b<=1; b++)
//                {
//                    TD[ a*2 + b ] = TDE[ a*4 + b*2 + 0 ] + TDE[ a*4 + b*2 + 1 ];
//                    TE[ a*2 + b ] = TDE[ a*4 + 0*2 + b ] + TDE[ a*4 + 1*2 + b ];                            
//                }
            
            //Calculate CT from CTE
            for(int a=0; a<=1; a++)
                for(int b=0; b<=1; b++)
                    CT[ a*2 + b ] = CTE[ a*4 + b*2 + 0 ] + CTE[ a*4 + b*2 + 1 ];
                            
                            
//            //Calculate T from TE
//            for(int a=0; a<=1 ;a++)   
//                T[a] = TE[ a*2 + 0 ] + TE[ a*2 + 1 ];
                            
            //Calculate C[] from CT
            for (int i=0; i<=1; i++)
                C[i] = CT[ i*2 + 0 ] + CT[ i*2 + 1 ];
             
//            //There is no count for T0
//            T[0] = 0.0; 
//            //There is no count for T0ge0, T0ge1
//            TE[0]=TE[1] = 0.0;
            //There is no count for C0T0, C1T0, ....C15T0
            for (int i=0; i<=1; i++)
                CT[ i*2 + 0 ] =  0.0;
            //There is no count for C0T0E0, C0T0E1, C1TOE0, C1T0E1,.....C15T0E0, C15T0E1
            for (int i=0; i<=1; i++)
                CTE[ i*4 + 0*2 + 0 ] = CTE[ i*4 + 0*2 + 1 ] = 0.0;
              //CTE[0] = CTE[1] = CTE[4] = CTE[5] = 0.0;
              

            
//           /**************************************************
//            * Calculate lnData NOT considering cancer type 
//            */
//            float lnData = 0.0;
//            
//            float TFscore ;
//            float DFscore ;
//            float pGT1GE1 ;
//            float pGT0GE1 ;
//     
//            if(curGTIndx == 0)
//            {
//                //TFscore = calcA0Fscore(T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0);
//                TFscore = calcA0Fscore(T[1],  TE[3], TE[2], T[0],  TE[1], TE[0]);
//            }
//            else 
//            {
//                //TFscore = calcFscore( T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0 );
//                TFscore = calcFscore( T[1],  TE[3], TE[2], T[0],  TE[1], TE[0] );
//            }
//
//            //float DFscore = calcFscore( D1, D1ge1, D1ge0, D0, D0ge1, D0ge0 );
//            DFscore = calcFscore( TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0] );
//
//            lnData = TFscore + DFscore + lntumorMutPriors[gt];
//
//            if(gt == 0)
//            {
//                //pGT1GE1 = (ALPHANULL + T1ge1) / (ALPHANULL + ALPHANULL + T1);
//                //pGT0GE1 = (ALPHANULL + D0ge1 + D1ge1) / (ALPHANULL + ALPHANULL + nTumors - T1);
//                pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]);
//                pGT0GE1 = (ALPHANULL + TDE[1] + TDE[3]) / (ALPHANULL + ALPHANULL + nTumors - T[1]);
//            }
//            else
//            {
//                //pGT1GE1 = (ALPHAIJK11 + T1ge1) / (ALPHAIJK11 + ALPHAIJK10 + T1);
//                //pGT0GE1 = (ALPHAIJK01 + D0ge1 + D1ge1) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T1);       
//                pGT1GE1 = (ALPHAIJK11 + TE[3]) / (ALPHAIJK11 + ALPHAIJK10 + T[1]);
//                pGT0GE1 = (ALPHAIJK01 + TDE[1] + TDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T[1]);                      
//            }
//
//            if(pGT1GE1 <= pGT0GE1)
//            {
//                lnData = -FLT_MAX;
//            }
//            
//            //save lnData 
//            tumorPosteriorMatrix[gt* nGE + ge] = lnData; 
//            
          
            /*****************************************************************
             *Calculate lnData for each cancel type , save the lnData into lnDataC[]
             */
 
  
            float lnDataC[2] = {0.0};
            
            float TFscoreC[2] = {0.0};
            float DFscoreC[2] = {0.0};
            float pGT1GE1C[2] = {0.0};
            float pGT0GE1C[2] = {0.0};             
            //float lnDataC1 = 0.0;

            for (int i=0; i<=1; i++)
            {
                if(curGTIndx == 0)
                {
                    //TFscore = calcA0Fscore(T1,   T1ge1,  T1ge0,   T0,   T0ge1,  T0ge0);//Original code
                    //corresponding to      C0T1,  C0T1E1, C0T1E0, C0T0,  C0T0E1, C0T0E0
                    //corresponding to      CT[1]] CTE[3]  CTE[2]  CT[0]] CTE[1]  CTE[0]]
                    //TFscore = calcA0Fscore(T[1],  TE[3], TE[2], T[0],  TE[1], TE[0]);//this corresponding to C=0 | GPU according to Original code
  
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
                //DFscore = calcFscore( CTD[5], CTDE[11], CTDE[10], CTD[4], CTDE[9], CTDE[8] );
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
                    lnDataC[i] = -FLT_MAX;
                }
                //save lnDataC[]]
                //tumorPosteriorMatrixC0[gt* nGE + ge] = lnDataC0;
                tumorPosteriorMatrixC[ (gt* nGE + ge ) + (nGT * nGE) * i ] = lnDataC[i];
                
            }
                
        }
        
//        float normalizer = 0.0;
        float normalizerC[2] = {0.0};
        
        for(int gt = 0; gt < nGT; gt++)
        {
            if(gt == 0)
            {
//                normalizer = tumorPosteriorMatrix[gt * nGE + ge];
                for (int i=0; i<=1; i++)
                    normalizerC[i] = tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i  ]; //each cancer type shifts nGT * nGE array index
                //normalizerC0 = tumorPosteriorMatrixC0[gt * nGE + ge]; 
                //normalizerC1 = tumorPosteriorMatrixC1[gt * nGE + ge];
            }
            else
            {
//                normalizer = logSum(normalizer, tumorPosteriorMatrix[gt * nGE + ge]);
                for (int i=0; i<=1; i++)
                    normalizerC[i] = logSum(normalizerC[i], tumorPosteriorMatrixC[ (gt * nGE + ge) + (nGT * nGE) * i ]);
                //normalizerC0 = logSum(normalizerC0, tumorPosteriorMatrixC0[gt * nGE + ge]);
                //normalizerC1 = logSum(normalizerC1, tumorPosteriorMatrixC1[gt * nGE + ge]);
            }
        }
        
        // finished populating a column of GTs with respect to a given GE, normalize so that the column sum to 1
        for (int gt = 0; gt < nGT; gt++)
        {
//            tumorPosteriorMatrix[gt * nGE + ge] = exp(tumorPosteriorMatrix[gt * nGE + ge] - normalizer);     
            for (int i=0; i<=1; i++)
               tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i] = exp(tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i] - normalizerC[i]);  
               
            //tumorPosteriorMatrixC0[gt * nGE + ge] = exp(tumorPosteriorMatrixC0[gt * nGE + ge] - normalizerC0);  
            //tumorPosteriorMatrixC1[gt * nGE + ge] = exp(tumorPosteriorMatrixC1[gt * nGE + ge] - normalizerC1); 
        }
        
 
//        float normT = 0.0; //without considering cancer type
        
        float normC[2] = {0.0}; //for each cancer type
        float normCT = 0.0; //combination of every cancer type normC[]
        
        

        
        for (int gt = 0; gt < nGT; gt++)
        {
            //normC0 = tumorPosteriorMatrixC0[gt * nGE + ge] ;
            //normC1 = tumorPosteriorMatrixC1[gt * nGE + ge];
//            normT =  tumorPosteriorMatrix[gt * nGE + ge];  
            //tumorPosteriorMatrix[gt * nGE + ge] = (normC0 * (nTumors - nC1)/nTumors + normC1 * nC1/nTumors + normT)/2;  
            normCT = 0.0;
            for (int i=0; i<=1; i++)
            {
                normC[i] = tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i];
                normCT += normC[i] * C[i]/nTumors;
            }          
                
//            tumorPosteriorMatrix[gt * nGE + ge] = (normCT + normT)/2;
            tumorPosteriorMatrix[gt * nGE + ge] = normCT;
                    
        }    
    }
}
  
/********** logSum *********************************************************/ 
/**
 * Evaluate Ln(x + y)
 * @param lnx ln(x)
 * @param lny ln(y)
 * @return ln(x + y)
 */
__device__ float logSum(float lnx, float lny){

    if(lnx == -INFINITY && lny == -INFINITY)
        return (lnx);

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
__device__ float calcFscore(float gt1,  float gt1ge1, float gt1ge0, 
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
__device__ float calcA0Fscore(float gt1,  float gt1ge1, float gt1ge0, 
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

