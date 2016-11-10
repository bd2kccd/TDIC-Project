 /*

 */

#include <math.h>
//#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <algorithm>
#include "TDIC_GeLoop_Cuda.h"
#include "PanCanTDIC_GeLoop_Cuda.h"

using namespace std;


 
/********** TDIC_GeLoop *********************************************************/ 
/**
 * TDIC at ge loop send to GPU for parallel computation
 * @param int d_nGE   Length of tumorGeIndieces. Note tumorGeIndices is a pointer now, not a vector, no way to use .size() to get nGE
 * @param int d_nGT   Length of tumorGtIndieces.
 * @param int nTumors   Number of total tumors
 * @param int *tumorGeIndices   Indices of Ge of given tumor
 * @param int *tumorGtIndices   Indices of Gt of given tumor
 * @param int *tumorGlobDriverIndices   Index of global driver gt for each ge of a given tumor
 * @param int d_maxDriverSize  maximal driver size
 * @param bool *gtDataMatrix    gt-by-tumor matrix represented as a 1-D array
 * @param bool *geDataMatrix    ge-by-tumor matrix represented as a 1-D array
 * @param float *lntumorMutPriors    Prior probability that an gt is a driver in a tumor
 * @param float *tumorPosteriorMatrix      nTumorGts-by-nTUmorGes matrix represented as a 1-D array. Store the posterior probability that a Gt causes a Ge. This param will be transfered out to CPU for outputting to files
 */   
__global__ void TDIC_GeLoop(int d_nGE, int d_nGT, int d_nTumors,  int* tumorGeIndices,  int* tumorGtIndices, int* tumorGlobDriverIndices, int d_maxDriverSize, bool* gtDataMatrix, bool* geDataMatrix, float* lntumorMutPriors, float* tumorPosteriorMatrix)
{
    
 float  ALPHANULL = 1.0;
 float  ALPHAIJK00 = 2.0;
 float  ALPHAIJK01 = 1.0;
 float  ALPHAIJK10 = 1.0;
 float  ALPHAIJK11 = 2.0;   
 
    int ge = threadIdx.x +  blockDim.x * blockIdx.x;//pass threadId to ge variable
     
    if (ge<d_nGE)
    {
 
        float normalizer = 0;
        int curGeIndx = tumorGeIndices[ge];
        int rowStartForGE = curGeIndx * d_nTumors; 
        
        // find the globDriver for this give ge   
//        int curGDriverIndx = tumorGlobDriverIndices[ge]; //curGDriverIndx is found by ge indx        
//        int rowStartForGlobDriver = curGDriverIndx * d_nTumors;
           
//        vector<int> curGDriverIndxes = tumorGlobDriverIndices[ge]; //curGDriverIndx is found by ge indxes        
//        vector<int> rowStartForGlobDrivers;
//        for(int i = 0; i < curGDriverIndxes.size(); i++)
//            rowStartForGlobDrivers.push_back(curGDriverIndxes[i]*nTumors);        
        
        // loop through each GT in the tumor
        for (int gt = 0; gt < d_nGT; gt++)
        {         
            
            float T[2] = {0.0};
            float TE[4] = {0.0};
            float TD[4] = {0.0};
            float TDE[8] = {0.0};
            
            
            int curGTIndx = tumorGtIndices[gt];

            int gtRowStart = curGTIndx * d_nTumors;
  
            int gdIndex, rowStartForGD;
            for(int t = 0; t < d_nTumors; t++)
            {
                int tVal = gtDataMatrix[gtRowStart + t];
                int eVal = geDataMatrix[rowStartForGE + t];
//                int dVal = gtDataMatrix[rowStartForGlobDriver + t\n
                //find dval
                int dVal = 0;
                for (int i = 0; i < d_maxDriverSize; i++)
                {
                    gdIndex = tumorGlobDriverIndices[ge*d_maxDriverSize + i];
                    rowStartForGD = gdIndex * d_nTumors;
                    dVal = dVal || gtDataMatrix[rowStartForGD + t];
                    if(dVal == 1)
                        break;
                }
                 
                TDE[tVal*4+dVal*2+eVal]++;
            }

            //TD[0xTD] T is the gt value, D is the global driver value
            //e.g. TD[2] means TD[ox10] save the count when T=1 D=0
            TD[0] = TDE[0] + TDE[1]; //T0D0 = T0D0E0 + T0D0E1 
            TD[1] = TDE[2] + TDE[3]; //T0D1 = T0D1E0 + T0D1E1
            TD[2] = TDE[4] + TDE[5]; //T1D0 = T1D0E0 + T1D0E1 
            TD[3] = TDE[6] + TDE[7]; //T0D1 = T1D1E0 + T1D1E1 
            //TE[0xTE]] T is the gt value, E is the ge value
            //e.g. TE[3] means TE[0x11] save the count when T=1 and E=1
            TE[0] = TDE[0] + TDE[2]; //T0E0 = T0D0E0 + T0D1E0
            TE[1] = TDE[1] + TDE[3]; //T0E1 = T0D0E1 + T0D1E1 
            TE[2] = TDE[4] + TDE[6]; //T1E0 = T1D0E0 + T1D1E0
            TE[3] = TDE[5] + TDE[7]; //T1E1 = T1D0E1 + T1D1E1
            //T[0xT] T is the gt value
            //e.g. T[1] save the count when gt value T = 1 
            T[0] = TE[0] + TE[1]; //T0 = T0E0 + T0E1
            T[1] = TE[2] + TE[3]; //T1 = T1E0 + T1E1  
    
            
            //Therr is no count for T0ge0, T0ge1 and T0
            TE[0]=TE[1] = 0.0;
            T[0] = 0.0;
                    
                    
            float TFscore;
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
            float DFscore = calcFscore( TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0] );

            float lnData = TFscore + DFscore + lntumorMutPriors[gt];
            

            tumorPosteriorMatrix[gt * d_nGE + ge] = lnData;
           

            float pGT1GE1, pGT0GE1;
            if(gt == 0)
            {
                //pGT1GE1 = (ALPHANULL + T1ge1) / (ALPHANULL + ALPHANULL + T1);
                //pGT0GE1 = (ALPHANULL + D0ge1 + D1ge1) / (ALPHANULL + ALPHANULL + d_nTumors - T1);
                pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]);
                pGT0GE1 = (ALPHANULL + TDE[1] + TDE[3]) / (ALPHANULL + ALPHANULL + d_nTumors - T[1]);
            }
            else
            {
                //pGT1GE1 = (ALPHAIJK11 + T1ge1) / (ALPHAIJK11 + ALPHAIJK10 + T1);
                //pGT0GE1 = (ALPHAIJK01 + D0ge1 + D1ge1) / (ALPHAIJK01 + ALPHAIJK00 + d_nTumors - T1);       
                pGT1GE1 = (ALPHAIJK11 + TE[3]) / (ALPHAIJK11 + ALPHAIJK10 + T[1]);
                pGT0GE1 = (ALPHAIJK01 + TDE[1] + TDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + d_nTumors - T[1]);                      
            }

   
            if(pGT1GE1 <= pGT0GE1)
            {
                tumorPosteriorMatrix[gt* d_nGE + ge] = -FLT_MAX;
            }
        }

        for(int gt = 0; gt < d_nGT; gt++)
        {
            if(gt == 0)
            {
                normalizer = tumorPosteriorMatrix[gt * d_nGE + ge];
            }
            else
            {
                normalizer = logSum(normalizer, tumorPosteriorMatrix[gt * d_nGE + ge]);
            }
        }
        
        // finished populating a column of GTs with respect to a given GE, normalize so that the column sum to 1
        for (int gt = 0; gt < d_nGT; gt++)
            tumorPosteriorMatrix[gt * d_nGE + ge] = exp(tumorPosteriorMatrix[gt * d_nGE + ge] - normalizer);  
        
       
    }
    
}

