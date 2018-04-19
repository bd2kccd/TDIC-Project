/* 
 * File:   TDIC.cpp
 * Author: Kevin Lu
 * 
 * Created on April 25, 2015, 6:13 PM
 */

using namespace std;
#include "TDIC.h"

#include <fstream>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <string>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <float.h>
#include <algorithm>


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

void TDIC(GTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        string>& mapGlobDrivers, vector<int>  tumorGtIndices, vector<double> gtPriorsTumor, vector<int> tumorGeIndices, const double v0, vector<double>& tumorPosteriorMatrix){
      
    
    //get Mat pointer
    bool* gtDataMatrix = gtMatrix.getMatPtr();
    bool* geDataMatrix = geMatrix.getMatPtr(); 
    
//    //calculate lntumorPriors
//    vector<double> lntumorMutPriors;
//    gtMatrix.calcLnTumorPriors(tumorGtIndices, v0, lntumorMutPriors);

    // Find the global drivers corresponding to the ge indices
    vector<int> tumorGlobDriverIndices;
    if (!getDEGGlobDriverIndices(gtMatrix, geMatrix, mapGlobDrivers, tumorGeIndices, tumorGlobDriverIndices))
    {
        cout << "Error occurred when retrieving global drivers";
    }
 

    //get nGT, nGE, nTumores
    unsigned int nGT = tumorGtIndices.size();
    unsigned int nGE = tumorGeIndices.size();
    int nTumors = gtMatrix.getNTumors();
    if (nTumors != geMatrix.getNTumors()) // number of tumors should agree between gtMatrix and geMatrix
    {
        cerr << "Bug: gtMatrix and geMatrix contain different number of tumors ";
    }

        
 
    // loop through each GE
    #pragma omp parallel for
    for(unsigned int ge = 0; ge < nGE; ge++)
    {
        double normalizer = 0;
        unsigned int curGeIndx = tumorGeIndices[ge];
        unsigned int rowStartForGE = curGeIndx * nTumors; 
        
        // find the globDriver for this give ge   
        unsigned int curGDriverIndx = tumorGlobDriverIndices[ge]; //curGDriverIndx is found by ge indx        
        unsigned int rowStartForGlobDriver = curGDriverIndx * nTumors;
        
        // loop through each GT in the tumor
        for (unsigned int gt = 0; gt < nGT; gt++)
        {
            // statistics associated with current T and global driver only 
            //float T1 = 0.0,   T1ge1 = 0.0, T1ge0 = 0.0, T0 = 0.0, T0ge1 = 0.0, T0ge0 = 0.0; 
            //float D1 = 0.0,   D1ge1 = 0.0, D1ge0 = 0.0, D0 = 0.0, D0ge1 = 0.0, D0ge0 = 0.0;
        
            float T[2] = {0.0};
            float TE[4] = {0.0};
            float TD[4] = {0.0};
            float TDE[8] = {0.0};
          
            int curGTIndx = tumorGtIndices[gt];
            int gtRowStart = curGTIndx * nTumors;
            
            double gtPrior;

            if (gtPriorsTumor[gt] == 0)
            {
                gtPrior = -FLT_MAX;
            }
            else
            {
                gtPrior = log(gtPriorsTumor[gt]);
            }
            
            for(int t = 0; t < nTumors; t++)
            {
                int tVal = gtDataMatrix[gtRowStart + t];
                int eVal = geDataMatrix[rowStartForGE + t];
                int dVal = gtDataMatrix[rowStartForGlobDriver + t];
               
//                // T, TE, TD can all be calculated from TDE
//                   T[tVal]++; 
//                   TE[tVal*2+eVal]++;
//                   TD[tVal*2+dVal]++;  
                

                //TDE[0xTDE] T is the gt value, D is the global driver value, E is the ge value
                //e.g. TDE[7] means TDE[0x110] save the count when T=1 and D=1 and E=1
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
              /*  
                //if GT = 1 at current tumor, update states for GT = 1 
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
            */
            
            //There is no count for T0ge0, T0ge1 and T0
            TE[0]=TE[1] = 0.0;
            T[0] = 0.0;
                    
                    
            double TFscore;
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
            double DFscore = calcFscore( TD[1], TDE[3], TDE[2], TD[0], TDE[1], TDE[0] );

            double lnData = TFscore + DFscore + gtPrior;

            tumorPosteriorMatrix[gt * nGE + ge] = lnData;

            double pGT1GE1, pGT0GE1;
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
        
        //for test only
//        float sumT = 0.0;
//        for (unsigned int gt = 0; gt < nGT; gt++)
//            sumT += tumorPosteriorMatrix[gt * nGE + ge] ;
//        if (fabs(sumT - 1) > 0.001)
//            cout << "sumT != 1 in Gene " << geNames[ge] << " difference is " << sumT - 1 << "\n";
        //test end
        
    }
}

/**
 * This function write TDIC or PanCanTDIC results to a csv file
 * @param tumorPosteriorMatrix: a vector contains the results of the TDIC which will be passed into this function to be saved to file
 * @param curTumorName: a string contains the current tumor name of the output
 * @param gtNames: a vector contains the gt gene names of the outputs
 * @param geNames: a vector contains the ge gene names of the outputs.
 * @param outPath: a string contains the file path of the output files
 */
void TDIC_Output(vector<double>& tumorPosteriorMatrix, string curTumorName, vector<string> gtNames, vector<string> geNames, string outPath){
    // save results to file
    string outFileName;
    if (*outPath.end() != '/')
    {
        outFileName = outPath + "/" +  curTumorName + ".csv";
    }
    else
    {
        outFileName = outPath + curTumorName + ".csv";
    }
        
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
    int nGE = geNames.size();
    for(int i = 0; i < nGE; i++)
    {
        outFile << "," << geNames[i];
    }
    outFile << "\n";
    
    int nGT = gtNames.size();
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
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end()); //added 4/14/16
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); //added 4/14/16
            fields = split(line, ',');
            globDriverMap.insert(std::pair<string, string>(fields.at(0), fields.at(1)));
                
        }
        inFileStream.close();
    }
    catch (ifstream::failure e) {
        cout << "Fail to open file " << fileName;
        return false;
    } 
   
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


/********** logSum *********************************************************/ 
/**
 * Evaluate Ln(x + y)
 * @param lnx ln(x)
 * @param lny ln(y)
 * @return ln(x + y)
 */
double logSum(double lnx, double lny){
    double maxExp = -4950.0;

    if(lny > lnx){                
        double tmp = lnx;
        lnx = lny;
        lny = tmp;
    }

    double lnyMinusLnX = lny - lnx;
    double lnXplusLnY;

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
double calcFscore(float gt1,  float gt1ge1, float gt1ge0, 
    float gt0, float gt0ge1, float gt0ge0 )
{
    // Calculation of Equation 7    
    double glnNi0 = lgamma(ALPHAIJK00 + ALPHAIJK01) - lgamma(gt0 + ALPHAIJK00 + ALPHAIJK01);
    double glnNi1 = lgamma(ALPHAIJK10 + ALPHAIJK11) - lgamma(gt1 + ALPHAIJK10 + ALPHAIJK11);

    double fscore = glnNi0 + glnNi1;   
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
double calcA0Fscore(float gt1,  float gt1ge1, float gt1ge0, 
    float gt0, float gt0ge1, float gt0ge0 )
{

    // Calculation of Equation 7    
    double glnNi0 = lgamma( ALPHANULL + ALPHANULL) - lgamma(gt0 + ALPHANULL + ALPHANULL);
    double glnNi1 = lgamma(ALPHAIJK10 + ALPHANULL) - lgamma(gt1 + ALPHANULL + ALPHANULL);

    double fscore = glnNi0 + glnNi1;   
    fscore += lgamma(gt0ge0 + ALPHANULL) - lgamma(ALPHANULL);
    fscore += lgamma(gt0ge1 + ALPHANULL) - lgamma(ALPHANULL);
    fscore += lgamma(gt1ge0 + ALPHANULL) - lgamma(ALPHANULL);
    fscore += lgamma(gt1ge1 + ALPHANULL) - lgamma(ALPHANULL);

    return (fscore);
}

/**
 * This function read in data input file, cancer type coding table file,  return SGA gene names, DEG gene names, and tumor name of each tumor
 * @param inputFileName(in):  A string contains the name of the customer tumor input file 
 *                       It is a txt file having 4 rows which shows as following. 
 *                       TumorID
 *                       CancerType
 *                       SGA genes list
 *                       DEG genes list
 * 
 * @geGetNames      a vector of a vector contains the input data gt gene names of each input tumor
 * @geGeneNames     a vector of a vector contains the input data ge gene names of each input tumor      
 * @TumorNames      a vector contains each tumor name  
 * @tumorCanTypes   a vector contains cancer type of each tumor
 * @cancerTypeTable a string which is the file name of the cancer type coding table for our TDIC. 
 */
void TDIC_Load(string inputFileName, vector< vector<string> >& gtGeneNames, vector< vector<double> >& gtGenePriors, vector< vector<string> >& geGeneNames, vector<string>& tumorNames, vector<int>& tumorCanTypes, string cancerTypeTable)
{
    std::stringstream ss;
    std::string line;
    ifstream inFileStream;   
    
    //Read in cancer type code table and saved in canTypeMap which is Cancer Type Name -> Code (1..16)
    try{
        inFileStream.open(cancerTypeTable.c_str());  
        if ( (inFileStream.rdstate() & std::ifstream::failbit ) != 0 )
        {
            std::cerr << "Error opening file" << cancerTypeTable <<" in TDIC_Load function, quit.\n";
            inFileStream.close();
            exit(EXIT_FAILURE);
        }
    }
    catch (...) { //std::ifstream::failure e
        cerr << "Fail to open file " << inputFileName;
        
    } 
   
    map<string, int> canTypeMap;
    while (getline(inFileStream, line)){
        if (!line.empty() && line[line.size() - 1] == '\r')
            line.erase(line.size() - 1);
            
        stringstream ss (line);
        string tmp1,tmp2;
        getline(ss, tmp1, ',');
        getline(ss, tmp2, ',');
        
        //remove the empty space of tmp1 and tmp2
        tmp1.erase(remove_if(tmp1.begin(), tmp1.end(), ::isspace),tmp1.end());
        tmp2.erase(remove_if(tmp2.begin(), tmp2.end(), ::isspace),tmp2.end());
        //tmp1 is cancer type, transform to upcase
        transform(tmp1.begin(), tmp1.end(), tmp1.begin(), ::toupper);
        
        canTypeMap[tmp1] = atoi(tmp2.c_str());
    }
    
    inFileStream.close();
   
    
    //Read in customer tumor input files
    try{
        inFileStream.open(inputFileName.c_str());  
        if ( (inFileStream.rdstate() & std::ifstream::failbit ) != 0 )
        {
            std::cerr << "Error opening file " << inputFileName << " in TDIC_Load function, quit.\n";
            inFileStream.close();
            exit(EXIT_FAILURE);
        }
    }
    catch (...) { //std::ifstream::failure e
        cerr << "Fail to open file " << inputFileName;
        
    } 
   
    
    int rowInFile = 0;
    vector<string> tumorGtNames; 
    vector<string> tumorGeNames;
    vector<double> tumorPriors;
    while (getline(inFileStream, line)){
        //remove '\r' in the line
        if (!line.empty() && line[line.size() - 1] == '\r')
            line.erase(line.size() - 1);
//        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        stringstream ss (line);
        string tmp;
        rowInFile ++;
        if (rowInFile%5 == 1){ //tumor name
            getline(ss, tmp, ',');
            tumorNames.push_back(tmp);
            tumorGtNames.clear(); 
            tumorGeNames.clear();
            tumorPriors.clear();
            tumorGtNames.push_back("A0");
        }

        else if (rowInFile%5 == 2) {//cancer type
            
//            getline(ss, tmp, ',');
//            int canType = canTypeMap[tmp.c_str()];
//            if (canType == 0){//if cancer type is not in our cancer type coding table, then map return 0
//                cerr << "Error: the name of cancer type " << tmp << " in the tumor input file is not correct./n";
//                exit(1);
//            }
//            tumorCanTypes.push_back( canType );
            
        }
        else if (rowInFile%5 == 3){ //Gt names
            while (getline(ss, tmp, ',')){
                tumorGtNames.push_back(tmp);
            }
        }
        else if (rowInFile%5 == 4){ //priors
            while (getline(ss, tmp, ',')){
                tumorPriors.push_back(atof(tmp.c_str()));
            }
        }
        else if (rowInFile%5 == 0){ //Ge names
            while (getline(ss, tmp, ',')){
                tumorGeNames.push_back(tmp);
            }
            
            //push tumorGtNames, tumorGeNames into gtGeneNames, geGeneNames;
            gtGeneNames.push_back(tumorGtNames);
            geGeneNames.push_back(tumorGeNames);
            gtGenePriors.push_back(tumorPriors);
        }

    }
    inFileStream.close();   

}


/**
 * This function pass in the gene names (gt/ge) of a tumor, and get their gene indices corresponding to TDIC ge/gt matrix 
 * @param inGeneNames(in): a vector contains the gene names of a input tumor
 * @param matrixGeneNames(in): a vector contains the gene Names of TDIC gt/ge matrix
 * @param gtIndices(out): a vector contains input genes indices corresponding to TDIC gtMatrix or geMatrix
 * 
 */
void getTumorGeneIndices(vector<string>& inGeneNames, vector<string>& matrixGeneNames, vector<int>& outGeneIndices)
{   
    for(int i = 0; i < inGeneNames.size(); i++)
    {
 
        int iFound = 0;
        for(int j = 0; j < matrixGeneNames.size(); j++)
        {

            if(inGeneNames[i].compare(matrixGeneNames[j]) == 0)
            {
                outGeneIndices.push_back(j);
                iFound = 1;
                break;
            }
        }
        if (iFound == 0)
            cout << "Customer input gene "<< inGeneNames[i] << " is not found in TDIC database\n";
    }
    
    if(inGeneNames.size() != outGeneIndices.size())
    {
        cout << "Some genes in the customer input file do not exist in TDIC database. \n";
    }
}