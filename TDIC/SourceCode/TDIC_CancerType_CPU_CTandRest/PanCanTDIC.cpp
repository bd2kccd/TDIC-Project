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
 * @param tumorID               The tumor to be processed
 */
void PanCanTDIC(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string,
        string>& mapGlobDrivers, const int tumorID, const string outPath, const double v0, int numCanType) {
    // initializations 

    string curTumorName = gtMatrix.getTumorNameById(tumorID);

    // Find the GTs and GEs of the current tumor
    vector<int> tumorGtIndices, tumorGeIndices;
    gtMatrix.findGeneWithOnesInTumor(tumorID, tumorGtIndices);
    geMatrix.findGeneWithOnesInTumor(tumorID, tumorGeIndices);

    bool* gtDataMatrix = gtMatrix.getMatPtr();
    bool* geDataMatrix = geMatrix.getMatPtr();

    // Find the global drivers corresponding to the
    vector<int> tumorGlobDriverIndices;
    //map <string, string> globalDriversMap;
    if (!getDEGGlobDriverIndices(gtMatrix, geMatrix, mapGlobDrivers, tumorGeIndices, tumorGlobDriverIndices)) {
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

    vector<double> lntumorMutPriors;
    gtMatrix.calcLnTumorPriors(tumorGtIndices, v0, lntumorMutPriors);

    double* tumorPosteriorMatrix = new double[nGT * nGE]();

    // loop through each GE
#pragma omp parallel for
    for (unsigned int ge = 0; ge < nGE; ge++) {
        int curGeIndx = tumorGeIndices[ge];
        int degRowStart = curGeIndx * nTumors;

        // find the globDriver for this given ge   
        int globalDriverIndx = tumorGlobDriverIndices[ge]; //curGDriverIndx is found by ge indx        
        int DriverRowStart = globalDriverIndx * nTumors;

        //Added according to python code
        int boolGlobalDriveAmongGT = 0; // flag to indicate if the global driver of current GE is in the GT limits
        int indxOfGlobalDriverInGTVec, indexOfMaxFscoreGTinGTVec;
        double maxFscore = -DBL_MAX;
        //        unsigned int curGeIndx = tumorGeIndices[ge];
        //        unsigned int rowStartForGE = curGeIndx * nTumors; 
        //        
        //        // find the globDriver for this given ge   
        //        unsigned int curGDriverIndx = tumorGlobDriverIndices[ge]; //curGDriverIndx is found by ge indx        
        //        unsigned int rowStartForGlobDriver = curGDriverIndx * nTumors;
        vector<int> canTypes = gtMatrix.getCanTypes();
        // loop through each GT in the tumor
        for (unsigned int gt = 0; gt < nGT; gt++) {
            int curGtIndx = tumorGtIndices[gt];
            int curGtRowStart = curGtIndx * nTumors;

            double fscore = calcPanCanFscore(gtDataMatrix, geDataMatrix, nTumors, lntumorMutPriors[gt], curGtRowStart, DriverRowStart,
                    degRowStart, curGtIndx, canTypes);

            tumorPosteriorMatrix[gt * nGE + ge] = fscore;


            //            // keep track the highest FScore and the score for the globa driver 
            //            if (fscore > maxFscore) {
            //                maxFscore = fscore;
            //                indexOfMaxFscoreGTinGTVec = gt;
            //            }
            //
            //            // check if the GlobalDriver for the GE is among the GTs
            //            if (curGtIndx == globalDriverIndx) {
            //                boolGlobalDriveAmongGT = 1;
            //                indxOfGlobalDriverInGTVec = gt;
            //            }
        } // end of looping through GTs

        //        // process the case of global drive is in the GT list
        //        if (boolGlobalDriveAmongGT){ 
        //            // switch the position of maxFscoreGT and global driver, recalculate fscore for the global driver in this tumor
        //            int curGtRowStart = DriverRowStart;
        //            int DriverRowStart = nTumors * tumorGtIndices[indexOfMaxFscoreGTinGTVec];
        //
        //            tumorPosteriorMatrix[nGE * indxOfGlobalDriverInGTVec + ge] = calcPanCanFscore(gtDataMatrix, geDataMatrix, nTumors,
        //                    lntumorMutPriors[indxOfGlobalDriverInGTVec], curGtRowStart, DriverRowStart,
        //                    degRowStart, globalDriverIndx, canTypes);
        //        }

        double normalizer = 0.0;

        for (int gt = 0; gt < nGT; gt++) {
            if (gt == 0) {
                normalizer = tumorPosteriorMatrix[gt * nGE + ge];
            } else {
                normalizer = logSum(normalizer, tumorPosteriorMatrix[gt * nGE + ge]);
            }
        }

        // finished populating a column of GTs with respect to a given GE, normalize so that the column sum to 1
        for (int gt = 0; gt < nGT; gt++) {
            tumorPosteriorMatrix[gt * nGE + ge] = exp(tumorPosteriorMatrix[gt * nGE + ge] - normalizer);
        }
    }

    // save results to file
    //vector<string> gtNames = gtMatrix.getGeneNames();

    string outFileName;
    if (*outPath.end() != '/') {
        outFileName = outPath + "/" + curTumorName + ".csv";
    } else {
        outFileName = outPath + curTumorName + ".csv";
    }


    //ofstream file;
    ofstream outFile;
    try {
        outFile.open(outFileName.c_str());
    } catch (ofstream::failure e) {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
    }

    //start writing CSV representation of TDIMatrix    
    //write column headers
    for (int i = 0; i < nGE; i++) {
        outFile << "," << geNames[i];
    }
    outFile << "\n";

    for (int i = 0; i < nGT; i++) {
        outFile << gtNames[i];
        for (int j = 0; j < nGE; j++) {
            outFile << "," << tumorPosteriorMatrix[i * nGE + j];
        }
        outFile << "\n";
    }
    outFile.close();

    delete [] tumorPosteriorMatrix;



}

/************ calcPanCanFscore ******************************************************************/
double calcPanCanFscore(const bool* gtDataMatrix, const bool* geDataMatrix, const int nTumors, const double curGtPrior,
        const int curGtRowStart, const int DriverRowStart, const int curGeRowStart, const int curGtIndx, vector<int> canTypes) {
    //Variables to save the count NOT considering the cancel type
    int numCanType = *max_element(canTypes.begin(), canTypes.end());
    double T[2] = {0.0};
    double TE[4] = {0.0};
    double TD[4] = {0.0};
    double TDE[8] = {0.0};

    //Variables to save the count considering the 16 cancel types

    double* C = new double[numCanType]();
    double* CT = new double[2 * numCanType]();
    double* CTE = new double[4 * numCanType]();
    double* CTD = new double[4 * numCanType]();
    double* CTDE = new double[8 * numCanType]();


    //count CTDE
    for (int t = 0; t < nTumors; t++) {

        int tumorCanType = canTypes[t];

        //int theSameCanType = (tumorCanType == compCanType);//if cancer type is the same, the value = 1 //this is for version cancer type the same or not

        int tVal = gtDataMatrix[curGtRowStart + t];
        int eVal = geDataMatrix[curGeRowStart + t];
        int dVal = gtDataMatrix[DriverRowStart + t];

        //                //Only need to count CTDE, all other variables can be calculated from CTDE

        //currentTumorCanType is from 1 to numCanType, the array starts from 0, so needs to minus 1
        CTDE[(tumorCanType - 1)*8 + tVal * 4 + dVal * 2 + eVal]++;
    }

    //Calculate CTE,CTD from CTDE
    for (int a = 0; a < numCanType; a++)
        for (int b = 0; b <= 1; b++)
            for (int c = 0; c <= 1; c++) {
                CTE[ a * 4 + b * 2 + c ] = CTDE[ a * 8 + b * 4 + 0 * 2 + c ] + CTDE[ a * 8 + b * 4 + 1 * 2 + c ];
                CTD[ a * 4 + b * 2 + c ] = CTDE[ a * 8 + b * 4 + c * 2 + 0 ] + CTDE[ a * 8 + b * 4 + c * 2 + 1 ];
            }
    //Calculate TDE from CTDE
    for (int a = 0; a <= 1; a++)
        for (int b = 0; b <= 1; b++)
            for (int c = 0; c <= 1; c++)
                for (int i = 0; i < numCanType; i++) //there are numCanType cancer types
                    TDE[ a * 4 + b * 2 + c ] += CTDE[ i * 8 + a * 4 + b * 2 + c ];

    //Calculate TD, TE from TDE
    for (int a = 0; a <= 1; a++)
        for (int b = 0; b <= 1; b++) {
            TD[ a * 2 + b ] = TDE[ a * 4 + b * 2 + 0 ] + TDE[ a * 4 + b * 2 + 1 ];
            TE[ a * 2 + b ] = TDE[ a * 4 + 0 * 2 + b ] + TDE[ a * 4 + 1 * 2 + b ];
        }

    //Calculate CT from CTE
    for (int a = 0; a < numCanType; a++)
        for (int b = 0; b <= 1; b++)
            CT[ a * 2 + b ] = CTE[ a * 4 + b * 2 + 0 ] + CTE[ a * 4 + b * 2 + 1 ];


    //Calculate T from TE
    for (int a = 0; a <= 1; a++)
        T[a] = TE[ a * 2 + 0 ] + TE[ a * 2 + 1 ];

    //Calculate C[] from CT
    for (int i = 0; i < numCanType; i++)
        C[i] = CT[ i * 2 + 0 ] + CT[ i * 2 + 1 ];

    //There is no count for T0
    T[0] = 0.0;
    //There is no count for T0ge0, T0ge1
    TE[0] = TE[1] = 0.0;
    //There is no count for CT0
    for (int i = 0; i < numCanType; i++)
        CT[ i * 2 + 0 ] = 0.0;
    //There is no count for CT0E0, CT0E1
    for (int i = 0; i < numCanType; i++)
        CTE[ i * 4 + 0 * 2 + 0 ] = CTE[ i * 4 + 0 * 2 + 1 ] = 0.0;
    // The following clearance for the storage of Crest count
    //There is no count for CT1D0, CT1D1
    for (int i = 0; i < numCanType; i++)
        CTD[ i * 4 + 1 * 2 + 0] = CTD[ i * 4 + 1 * 2 + 1] = 0.0;
    //There is no count for CT1D0E0, CT1D0E1,CT1D1E0, CT1D1E1
    for (int i = 0; i < numCanType; i++)
        for (int a = 0; a <= 1; a++)
            for (int b = 0; b <= 1; b++)
                CTDE[ i * 8 + 1 * 4 + a * 2 + b ] = 0.0;



    //    for (int c = 0; c < numCanType; c++)
    //        printf("gtIndex=%d, degRowStart=%d - %f,%f,%f\n",curGtIndx,curGeRowStart, CT[ c * 2 + 1 ],CTE[ 4 * c + 2], CTE[ 4 * c + 3]);
    /**************************************************
     * Calculate lnData NOT considering cancer type 
     */
    double lnData = 0.0;

    double TFscore;
    double DFscore;
    double pGT1GE1;
    double pGT0GE1;

    // printf("AAA %f,%f,%f,%d\n", TE[3], T[1], TDE[1]+TDE[3], nTumors);
    if (curGtIndx == 0) {
        //pGT1GE1 = (ALPHANULL + T1ge1) / (ALPHANULL + ALPHANULL + T1);
        //pGT0GE1 = (ALPHANULL + D0ge1 + D1ge1) / (ALPHANULL + ALPHANULL + nTumors - T1);
        pGT1GE1 = (ALPHANULL + TE[3]) / (ALPHANULL + ALPHANULL + T[1]);
        pGT0GE1 = (ALPHANULL + TDE[1] + TDE[3]) / (ALPHANULL + ALPHANULL + nTumors - T[1]);
    } else {
        //pGT1GE1 = (ALPHAIJK11 + T1ge1) / (ALPHAIJK11 + ALPHAIJK10 + T1);
        //pGT0GE1 = (ALPHAIJK01 + D0ge1 + D1ge1) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T1);       
        pGT1GE1 = (ALPHAIJK11 + TE[3]) / (ALPHAIJK11 + ALPHAIJK10 + T[1]);
        pGT0GE1 = (ALPHAIJK01 + TDE[1] + TDE[3]) / (ALPHAIJK01 + ALPHAIJK00 + nTumors - T[1]);
    }

    if (pGT1GE1 <= pGT0GE1 and curGtIndx != 0) {
        lnData = -DBL_MAX;
        delete[] C;
        delete[] CT;
        delete[] CTE;
        delete[] CTD;
        delete[] CTDE;
        
        return lnData;
    }


    if (curGtIndx == 0) {
        //TFscore = calcA0Fscore(T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0);
        TFscore = calcA0Fscore(T[1], TE[3], TE[2], T[0], TE[1], TE[0]);
    } else {
        //TFscore = calcFscore( T1,  T1ge1, T1ge0, T0,  T0ge1, T0ge0 );
        TFscore = calcFscore(T[1], TE[3], TE[2], T[0], TE[1], TE[0]);
    }

    //double DFscore = calcFscore( D1, D1ge1, D1ge0, D0, D0ge1, D0ge0 );
    DFscore = calcFscore(TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0]);

    lnData = TFscore + DFscore + curGtPrior + log(0.5); //TCIRATIO=0.5

    /*****************************************************************
     *Calculate lnData considering each cancel type , save the lnData into lnDataC[]
     */
    double lnDataC = 0.0;
    double TFscoreC = 0.0;
    double DFscoreC = 0.0;
    for (int c = 0; c < numCanType; c++) 
    {    //Calculate the count for the rest of C
        for (int j = 0; j < numCanType; j++) {
            if (j != c) {
                //calculate CT0,CT0E1, CT0E0 count summing rest of canType of CT1,CT1E1, CT1E0 
                CT[ 2 * c + 0 ] += CT[ 2 * j + 1 ];
                CTE[ 4 * c + 2 * 0 + 1] += CTE[ 4 * j + 2 * 1 + 1];
                CTE[ 4 * c + 2 * 0 + 0] += CTE[ 4 * j + 2 * 1 + 0];
                //summing rest of CT0D1, CT0D0(T=0) of canType and saved into T1D1, T1D0(T=1)
                CTD[ 4 * c + 2 * 1 + 0] += CTD[ 4 * j + 2 * 0 + 0];
                CTD[ 4 * c + 2 * 1 + 1] += CTD[ 4 * j + 2 * 0 + 1];
                //summing rest of CT0DE(T=0) of canType and saved into CT1DE(T=1)
                for (int a = 0; a <= 1; a++)
                    for (int b = 0; b <= 1; b++)
                        CTDE[ 8 * c + 4 * 1 + 2 * a + b ] += CTDE[ 8 * j + 4 * 0 + 2 * a + b ];
            }
        }

//        //calculate Fscore for cancer type C, summing T0, T0E of rest of canType as T1, T1E
//        if (curGtIndx == 0) {
//            //TFscore = calcA0Fscore(T1,   T1ge1,  T1ge0,   T0,   T0ge1,  T0ge0);//Original code
//            //corresponding to      C0T1,  C0T1E1, C0T1E0, C0T0,  C0T0E1, C0T0E0
//            //corresponding to      CT[1]] CTE[3]  CTE[2]  CT[0]] CTE[1]  CTE[0]]
//            //TFscore = calcA0Fscore(T[1],  TE[3], TE[2], T[0],  TE[1], TE[0]);//this corresponding to C=0 | GPU according to Original code
//            TFscoreC = calcA0Fscore(CT[ 2 * c + 1 ], CTE[ 4 * c + 3], CTE[ 4 * c + 2], CT[ 2 * c + 0 ], CTE[ 4 * c + 1], CTE[ 4 * c + 0]);
//        } else {
//            TFscoreC = calcFscore(CT[ 2 * c + 1 ], CTE[ 4 * c + 3], CTE[ 4 * c + 2], CT[ 2 * c + 0 ], CTE[ 4 * c + 1], CTE[ 4 * c + 0]);
//        }

        //calculate FscoreC for cancer type C, 
        if (curGtIndx == 0) {
            //TFscore = calcA0Fscore(T1,   T1ge1,  T1ge0,   T0,   T0ge1,  T0ge0);//Original code
            //corresponding to      C0T1,  C0T1E1, C0T1E0, C0T0,  C0T0E1, C0T0E0
            //corresponding to      CT[1]] CTE[3]  CTE[2]  CT[0]] CTE[1]  CTE[0]]
            //TFscore = calcA0Fscore(T[1],  TE[3], TE[2], T[0],  TE[1], TE[0]);//this corresponding to C=0 | GPU according to Original code
            TFscoreC = calcA0Fscore(CT[ 2 * c + 1 ], CTE[ 4 * c + 3], CTE[ 4 * c + 2],0,0,0);
        } else {
            TFscoreC = calcFscore(CT[ 2 * c + 1 ], CTE[ 4 * c + 3], CTE[ 4 * c + 2], 0,0,0);
        }

        //calculate FscoreCrest for rest of canType c, added to the FsscoreC
        if (curGtIndx == 0) {
            TFscoreC += calcA0Fscore( CT[ 2 * c + 0 ], CTE[ 4 * c + 1], CTE[ 4 * c + 0],0,0,0);
        } else {
            TFscoreC += calcFscore(CT[ 2 * c + 0 ], CTE[ 4 * c + 1], CTE[ 4 * c + 0],0,0,0);
        }


        /*
         * calculate  DfscoreC for T=0  
         */
        //float DFscore = calcFscore( D1,     D1ge1,    D1ge0,     D0,     D0ge1,     D0ge0 );//Original code
        //corresponding to           CT0D1  CT0D1E1  CT0D1E0   CT0D0   CT0D0E1   CT0D0E0
        //corresponding to           CTD[1]  CTDE[3] CTDE[2]   CTD[0]  CTDE[1]   CTDE[0]
        //DFscore = calcFscore( TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0] );//GPU according to Original code
        DFscoreC = calcFscore(CTD[ 4 * c + 1], CTDE[ 8 * c + 3], CTDE[8 * c + 2], CTD[ 4 * c + 0 ], CTDE[ 8 * c + 1 ], CTDE[ 8 * c + 0]);
        /*
         * calculate C0 DfscoreCrest for T=1 
         */
        //float DFscore = calcFscore( D1,     D1ge1,    D1ge0,     D0,     D0ge1,     D0ge0 );//Original code
        //corresponding to           CT1D1  CT1D1E1  CT1D1E0   CT1D0   CT1D0E1   CT1D0E0
        //corresponding to           CTD[3]  CTDE[7] CTDE[6]   CTD[2]  CTDE[5]   CTDE[4]
        DFscoreC += calcFscore(CTD[ 4 * c + 3], CTDE[ 8 * c + 7], CTDE[8 * c + 6], CTD[ 4 * c + 2 ], CTDE[ 8 * c + 5 ], CTDE[ 8 * c + 4]);
        
        lnDataC = TFscoreC + DFscoreC + curGtPrior + log((0.5) * C[c] / nTumors);

        lnData = logSum(lnData, lnDataC);
    }

    //delete dynamic arrays
    delete[] C;
    delete[] CT;
    delete[] CTE;
    delete[] CTD;
    delete[] CTDE;
    return lnData;
}





