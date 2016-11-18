/* 
 * File:   PanCanTDIC_Global.cpp
 *
 * 
 * 
 */

//#include "TDIC.h"
#include "GenePanCanGlobDriver.h"
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
 * @param outFileName       A string that contains the output file name include the path 
 * @param v0   A constant float
 */
//void PanCanTDIC(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
//        string>& mapGlobDrivers, const int tumorID, const string outPath, const float v0){
void GenePanCanGlobDriver(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix,  const string outFileName, const float v0){    
    
    

    int* gtDataMatrix = gtMatrix.getMatPtr();
    int* geDataMatrix = geMatrix.getMatPtr();

    const int nCanType = gtMatrix.getNumCanType();
    
    // Allocate memory for nGT x nGE matrix

    unsigned int nGT = gtMatrix.getNGenes();
    unsigned int nGE = geMatrix.getNGenes();
    
    int nTumors = gtMatrix.getNTumors();
    if (nTumors != geMatrix.getNTumors()) // number of tumors should agree between gtMatrix and geMatrix
    {
        cerr << "Bug: gtMatrix and geMatrix contain different number of tumors ";
    }


    //Get names of  GTs and GEs 
     
    vector<string> gtNames = gtMatrix.getGeneNames();
    vector<string> geNames = geMatrix.getGeneNames();
    
    //define lnGlobalPriors with size of nGT and initialized to 0
    vector<float> lnGlobalPriors (nGT,0); 
    gtMatrix.calcGlobalAlphaPriors(v0, lnGlobalPriors);

    float* tumorPosteriorMatrix = new float[nGT * nGE]();
    //nGT * nGE * nCanType will be over 20G, cause segmentation fault. Create 2 dimensional array instead
    //float* tumorPosteriorMatrixC = new float[nGT * nGE * nCanType](); 
    float **tumorPosteriorMatrixC = new float*[nCanType];
    for( int i=0; i<nCanType; i++ ) {
        tumorPosteriorMatrixC[i] = new float[nGT * nGE]();
//        for( int j=0; j<nGT * nGE; j++ ) {
//            tumorPosteriorMatrixC[i][j] = 0.0;
//        }
    }
 
  
    // loop through each GE
    #pragma omp parallel for
    for(unsigned int ge = 0; ge < nGE; ge++)
    {
        if (ge % 1 == 0)
            printf("TDIC processing %d th gene.\n", ge);
        //C[] saves the count numbers of the each cancer type
        float* C = new float[nCanType]();
   
        unsigned int rowStartForGE = ge * nTumors; 

        // loop through each GT in the tumor
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            
            //Variables to save the count NOT considering the cancel type
            float T[2] = {0.0};
            float TE[4] = {0.0};

            //Variables to save the count considering the nCanType cancel types
            float* CT = new float[2*nCanType]();
            float* CTE = new float[4*nCanType]();

            int gtRowStart = gt * nTumors;
            
            //count CTE
            for(int t = 0; t < nTumors; t++)
            {
                
                int tumorCanType = gtMatrix.getCanTypeByTumorId(t);
                
                int tVal = gtDataMatrix[gtRowStart + t];
                int eVal = geDataMatrix[rowStartForGE + t];
                
//                //Only need to count CTE, all other variables can be calculated from CTE

                //currentTumorCanType is from 1 to ..., the array starts from 0, so needs to minus 1
                CTE[(tumorCanType-1)*4+tVal*2+eVal]++;
            }
            

            //Calculate TE from CTE
            //combine nCantype of CTE into TE which has a*b combinations
            for(int a=0; a<=1; a++)
                for(int b=0; b<=1; b++)
                    for(int i=0; i<nCanType; i++)
                        TE[ a*2 + b ] += CTE[ i*4 + a*2 + b ];
                
                
            //Calculate CT from CTE
            //combine CTE which has E=0 and E=1 into CT which has a*B combinations
            for(int a=0; a<nCanType; a++)
                for(int b=0; b<=1; b++)
                    CT[ a*2 + b ] = CTE[ a*4 + b*2 + 0 ] + CTE[ a*4 + b*2 + 1 ];
                            
                            
            //Calculate T from TE
            for(int a=0; a<=1 ;a++)   
                T[a] = TE[ a*2 + 0 ] + TE[ a*2 + 1 ];
                            
            //Calculate C from CT
            for (int i=0; i<nCanType; i++)
                C[i] = CT[ i*2 + 0 ] + CT[ i*2 + 1 ];
             
//            //There is no count for T0
//            T[0] = 0.0; 
//            //There is no count for T0ge0, T0ge1
//            TE[0]=TE[1] = 0.0;
//            //There is no count for C0T0, C1T0, ....C15T0
//            for (int i=0; i<nCanType; i++)
//                CT[ i*2 + 0 ] =  0.0;
//            //There is no count for C0T0E0, C0T0E1, C1TOE0, C1T0E1,.....C15T0E0, C15T0E1
//            for (int i=0; i<nCanType; i++)
//                CTE[ i*4 + 0*2 + 0 ] = CTE[ i*4 + 0*2 + 1 ] = 0.0;
//              //CTE[0] = CTE[1] = CTE[4] = CTE[5] = 0.0;
              

            
           /**************************************************
            * Calculate lnData NOT considering cancer type 
            */
            float lnData = 0.0;
            
            float TFscore ;
            float DFscore ;
            float pGT1GE1 ;
            float pGT0GE1 ;
     
            if(gt == 0)
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
//            DFscore = calcFscore( TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0] );

//            lnData = TFscore + DFscore + lntumorMutPriors[gt];
            lnData = TFscore + lnGlobalPriors[gt];

            if(gt == 0)
            {
                //pGT1GE1 = (ALPHANULL + T1ge1) / (ALPHANULL + ALPHANULL + T1);
                //pGT0GE1 = (ALPHANULL + D0ge1 + D1ge1) / (ALPHANULL + ALPHANULL + nTumors - T1);
                pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]);
//                pGT0GE1 = (ALPHANULL + TDE[1] + TDE[3]) / (ALPHANULL + ALPHANULL + nTumors - T[1]);
                pGT0GE1 = (ALPHANULL + TE[1]) / (ALPHANULL + ALPHANULL + T[0]);
            }
            else
            {
                //pGT1GE1 = (ALPHAIJK11 + T1ge1) / (ALPHAIJK11 + ALPHAIJK10 + T1);
                //pGT0GE1 = (ALPHAIJK01 + D0ge1 + D1ge1) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T1);       
                pGT1GE1 = (ALPHAIJK11 + TE[3]) / (ALPHAIJK11 + ALPHAIJK10 + T[1]);
//                pGT0GE1 = (ALPHAIJK01 + TDE[1] + TDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T[1]);  
                pGT0GE1 = (ALPHAIJK01 + TE[1]) / (ALPHAIJK01 + ALPHAIJK00 + T[0]); 
            }

            if(pGT1GE1 <= pGT0GE1)
            {
                lnData = -FLT_MAX;
            }
            
            //save lnData 
            tumorPosteriorMatrix[gt* nGE + ge] = lnData; 
            
          
            /*****************************************************************
             *Calculate lnData for nCanType cancel types , save the lnData into lnDataC[]
             */
 
            float* lnDataC = new float[nCanType]();
            
            float* TFscoreC = new float[nCanType]();
//            float* DFscoreC = new float[nCanType]();
            float* pGT1GE1C = new float[nCanType]();
            float* pGT0GE1C = new float[nCanType]();             
            //float lnDataC1 = 0.0;

            for (int i=0; i<nCanType; i++)
            {
                if(gt == 0)
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
//                DFscoreC[i] = calcFscore( CTD[1+4*i], CTDE[3+8*i], CTDE[2+8*i], CTD[0+4*i], CTDE[1+8*i], CTDE[0+8*i] );

//                lnDataC[i]= TFscoreC[i] + DFscoreC[i] + lntumorMutPriors[gt];
                lnDataC[i] = TFscoreC[i] + lnGlobalPriors[gt];
                

                if(gt == 0)
                {
                    //pGT1GE1 = (ALPHANULL + T1ge1) / (ALPHANULL + ALPHANULL + T1);  //Original code
                    //pGT0GE1 = (ALPHANULL + D0ge1 + D1ge1) / (ALPHANULL + ALPHANULL + nTumors - T1); //Original code

                    //pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]); //GPU according to Original code
                    //pGT0GE1 = (ALPHANULL + TDE[1] + TDE[3]) / (ALPHANULL + ALPHANULL + nTumors - T[1]);//GPU according to Original code
                    //TDE[1] + TDE[3] = TE[1]; nTumores - T[1] = T[0]
                    //pGT1GE1 = (ALPHANULL + CTE[7]) / (ALPHANULL + ALPHANULL + CT[3]);//GPU for cancer type 1 larger
                    //pGT0GE1 = (ALPHANULL + CTDE[9] + CTDE[11]) / (ALPHANULL + ALPHANULL + nTumors - CT[3]);//GPU for cancer type 1 larger
                    pGT1GE1C[i] = (ALPHANULL + CTE[3+4*i]) / (ALPHANULL + ALPHANULL + CT[1+2*i]);
                    pGT0GE1C[i] = (ALPHANULL + CTE[1+4*i]) / (ALPHANULL + ALPHANULL + CT[0+2*i]);
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
                    pGT0GE1C[i] = (ALPHAIJK01 + CTE[1+4*i]) / (ALPHAIJK01 + ALPHAIJK00 + CT[0+2*i]);    
                }

                if(pGT1GE1C[i] <= pGT0GE1C[i])
                {
                    //tumorPosteriorMatrix[gt* nGE + ge] = -FLT_MAX;      Changed for cancel type calculation. Save lnData at the last step
                    lnDataC[i] = -FLT_MAX;
                }
                //save lnDataC[]]
//                tumorPosteriorMatrixC[ (gt* nGE + ge ) + (nGT * nGE) * i ] = lnDataC[i];
                tumorPosteriorMatrixC[i][ gt* nGE + ge  ] = lnDataC[i];
                
            }
        
            //delete dynamic arrays
            delete[] CT;
            delete[] CTE;
            delete[] lnDataC;
            delete[] TFscoreC;
            delete[] pGT1GE1C;
            delete[] pGT0GE1C;
        
        }
        
        float normalizer = 0.0;
        float* normalizerC = new float[nCanType]();
        
        for(unsigned int gt = 0; gt < nGT; gt++)
        {
            if(gt == 0)
            {
                normalizer = tumorPosteriorMatrix[gt * nGE + ge];
                for (int i=0; i<nCanType; i++)
//                    normalizerC[i] = tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i  ]; //each cancer type shifts nGT * nGE array 
                    normalizerC[i] = tumorPosteriorMatrixC[i][ gt * nGE + ge ]; //each cancer type shifts nGT * nGE array index
  
            }
            else
            {
                normalizer = logSum(normalizer, tumorPosteriorMatrix[gt * nGE + ge]);
                for (int i=0; i<nCanType; i++)
//                    normalizerC[i] = logSum(normalizerC[i], tumorPosteriorMatrixC[ (gt * nGE + ge) + (nGT * nGE) * i ]);
                    normalizerC[i] = logSum(normalizerC[i], tumorPosteriorMatrixC[i][ gt * nGE + ge ]);
            }
        }
        
        // finished populating a column of GTs with respect to a given GE, normalize so that the column sum to 1
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            tumorPosteriorMatrix[gt * nGE + ge] = exp(tumorPosteriorMatrix[gt * nGE + ge] - normalizer);     
            for (int i=0; i<nCanType; i++)
//               tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i] = exp(tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i] - normalizerC[i]);  
               tumorPosteriorMatrixC[i][ gt * nGE + ge ] = exp(tumorPosteriorMatrixC[i][ gt * nGE + ge ] - normalizerC[i]);  
               
        }
        
 
        float normT = 0.0; //without considering cancer type
        
        float* normC = new float[nCanType](); //for nCanType cancer type
        float normCT = 0.0; //combination of nCanType cancer type normC[]
        
        
        float sumT = 0.0; //for test purpose
        float sumCT = 0.0; //for test purpose
        
        
        for (unsigned int gt = 0; gt < nGT; gt++)
        {

            normT =  tumorPosteriorMatrix[gt * nGE + ge];  
            //tumorPosteriorMatrix[gt * nGE + ge] = (normC0 * (nTumors - nC1)/nTumors + normC1 * nC1/nTumors + normT)/2;  
            normCT = 0.0;
            for (int i=0; i<nCanType; i++)
            {
//                normC[i] = tumorPosteriorMatrixC[(gt * nGE + ge) + (nGT * nGE) * i];
                normC[i] = tumorPosteriorMatrixC[i][ gt * nGE + ge ];
                normCT += normC[i] * C[i]/nTumors;
            }          
                
            tumorPosteriorMatrix[gt * nGE + ge] = (normCT + normT)/2;
                    
//            sumCT += normCT;//for test purpose only
//            sumT += normT; //for test purpose only
        }    
//        if (fabs(sumCT-1) > 0.01) //for test purpose
//            cout << "bug: in PanCanTDIC.Cpp/PanCanTDIC sumC != 1 difference is " << sumCT-1 << "\n"; //for test purpose
//        if (fabs(sumT - 1) > 0.01) //for test purpose
//            cout << "bug: in PanCanTDIC.Cpp/PanCanTDIC sumT != 1 difference is " << sumT-1 << "\n"; //for test purpose
        
        //delete dynamic arrays
        delete[] C;
        delete[] normalizerC;
        delete[] normC;
    }
  

    // save results to file
 
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
    //delete tumorPosteriorMatrixC array
    for( int i=0; i<nCanType; i++ ) {
        delete[] tumorPosteriorMatrixC[i];
    }
    delete [] tumorPosteriorMatrixC;

}
