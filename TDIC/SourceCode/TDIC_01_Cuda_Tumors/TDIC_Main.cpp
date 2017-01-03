#include <cstdlib>
#include <unistd.h>
#include <getopt.h>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <fstream>
#include <omp.h>
#include <time.h>
#include <algorithm>
//#include "TDIC.h"
#include "PanCanTDIC.h"
//#include "TDIMatrix.h"
//#include "GTMatrix.h"
#include "PanCanGTMatrix.h"


using namespace std;


int main(int argc, char** argv) {
    
    float v0 = 0.1;
    int hasOpt;
    string gtFilePath, gtcFilePath, globalDriverPath, degFilePath, outPath, inputDataFile, cancerTypeTable ;
    int rowStart = -1;
    int rowEnd = -1;

    time_t t_start,t_end;
    time (&t_start);

    while((hasOpt = getopt(argc, argv, "hi:f:d:g:o:c:s:e:t:?")) != -1)
    {
        switch(hasOpt)
        {
            case 't':
                cancerTypeTable = optarg;
                break;
                
            case 'i':
                inputDataFile = optarg;
                break;

            case 'f':
                gtFilePath = optarg;
                break;
                
            case 'c':
                gtcFilePath = optarg;
                break;
   
            case 'd':
                degFilePath = optarg;
                break;

            case 'g':
                globalDriverPath = optarg ;
                break;
            
            case 'o':
                outPath = optarg;
                break;
            case 's':
                if(atoi(optarg) >= 0)       
                    rowStart = atoi(optarg);
                else
                {
                    cout << "rowStart given is less than zero. Exiting out.\n";
                    exit(1);
                }
                break;
            
            case 'e':
                rowEnd = atoi(optarg);
                if(rowEnd < 0)
                {
                    cout << "rowEnd given is less than zero. Exiting out.\n";
                    exit(1);                    
                }
                break;

            case 'h':
                cerr << "Usage: TDIC -i inputData -c inputGtMatrix -d inputGeMatrix -g inputGlobDriverDictionary -o pathForOutputResults \n";
                exit(1);
                break;
                
            case '?':
                if(optopt == 'i')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;
                }
                else if(optopt == 'g')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;          
                }
                else if(optopt == 'd')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;          
                }
  
                else if(isprint(optopt))
                {
                  cout << "Unknown option -" << optopt << ".\n";
                  return 0;
                }
                else
                {
                  cout << "Unknown option character.\n";
                  return 0;
                }

            default:
                cerr << "Usage: TDIC -i inputDataFile -c inputPancanGtMatrix [-f inputGtMatrix] -d inputGeMatrix -g inputGlobDriverDictionary -o pathForOutputResults -t cancerTypeCodingTable\n";
                abort();
        }
    }
     
    
    GTMatrix* gtMatrix;
    PanCanGTMatrix* panCanGtMatrix;

    if (!gtFilePath.empty() && !gtcFilePath.empty() )//input both GTMatrix and PanCanGTMatrix
    {
        //gtFileePath and gtcFilePath can not both exist, either process GTMatrix or PanCanGTMatrix
        cerr << "Can not input both GtMatrix and PanCanGtMatrix\n";
        exit(1);            
    }
    else if (gtFilePath.empty() && gtcFilePath.empty() )//both gtFileePath and gtcFilePath not exist
     {
        
        cerr << "Must input GtMatrix or PanCanGtMatrix \n";
        exit(1);            
    }    
    else if (inputDataFile.empty() )
    {
        cerr << "Must have input data file\n";
        exit(1);
    }
    else if (globalDriverPath.empty())
    {
        cerr << "Must input global driver file\n";
        exit(1);
    }
    else if (!gtFilePath.empty()) //input GTMatrix
    {          
        cerr << "Reading GT matrix: " << gtFilePath << "\n";
        gtMatrix = new GTMatrix(gtFilePath);
    }
    else //input PanCanGTMatrix
    {   cout << "Reading PanCanGT matrix: " << gtcFilePath << "\n";
        panCanGtMatrix = new PanCanGTMatrix(gtcFilePath);
    }   
     
    
    //read in GE matrices
    cout << "Reading GE matrix. " << degFilePath << "\n";
    TDIMatrix* geMatrix = new TDIMatrix(degFilePath);
    
    int nTumors = geMatrix->getNTumors();
    
    
    cout << "Reading global driver file.\n";
    map<string, string> globalDriverMap;
    parseGlobDriverDict(globalDriverPath, globalDriverMap);

    
    /*
     * read in input data file and loop through tumors to get tumor gt ge gene indices then run TDI
     */
   
    //Call TDIC_Load to get gt, ge names
    vector< vector <string> > gtGeneNames, geGeneNames;
    vector <string> tumorNames;
    vector<int> tumorCanTypes;
    TDIC_Load(inputDataFile, gtGeneNames, geGeneNames, tumorNames, tumorCanTypes, cancerTypeTable);

    int inputNtumors = tumorNames.size();
    if(rowStart < 0)
        rowStart = 0;
    if(rowEnd < 1 || rowEnd > inputNtumors )
        rowEnd = inputNtumors;

    if(rowStart >= rowEnd)
    {
        cout << "Given rowEnd index must bigger than given rowStart. Exiting out.\n";
        exit(1);
    }

    //Get TDIC gt/ge gene names
    vector<string> gtMatrixGeneNames, geMatrixGeneNames;
    if (!gtcFilePath.empty())
        gtMatrixGeneNames = panCanGtMatrix->getGeneNames();
    else
        gtMatrixGeneNames = gtMatrix->getGeneNames();
    geMatrixGeneNames = geMatrix->getGeneNames();
    
    
    //prepare non-tumor specific device variable for GPU kernel invoke
    //define the host variables
    vector<int> cancerTypes;
    bool *gtDataMatrix, *geDataMatrix;
    int numCanTypes;
    if (!gtcFilePath.empty()){
        gtDataMatrix = panCanGtMatrix->getMatPtr();
        cancerTypes = panCanGtMatrix->canTypes;
         //Get number of cancer types. Use max value of the cancer type as the number of cancer types
        numCanTypes = *max_element(cancerTypes.begin(), cancerTypes.end());
        cout <<  "There are total of " << numCanTypes << " cancer types.\n";
    }
    else 
        gtDataMatrix = gtMatrix->getMatPtr();
    
    geDataMatrix  = geMatrix->getMatPtr();
    
    //define and malloc device variables
    int  *d_cancerTypes;
    if (!gtcFilePath.empty()){
        cudaMalloc( (int**)&d_cancerTypes, nTumors*sizeof(int) );
    }
    
    bool *d_gtDataMatrix, *d_geDataMatrix;
    int numColT,numRowT;
    if (!gtcFilePath.empty()){
        numColT = panCanGtMatrix->nCol;
        numRowT = panCanGtMatrix->nRow;
    }
    else{
        numColT = gtMatrix->nCol;
        numRowT = gtMatrix->nRow;
    }
    cudaMalloc( (int**)&d_gtDataMatrix,numColT*numRowT*sizeof(bool) );

    int numColE = geMatrix->nCol;
    int numRowE = geMatrix->nRow;
    cudaMalloc( (int**)&d_geDataMatrix,numColE*numRowE*sizeof(bool) );

    //transfer data from host to device
    if (!gtcFilePath.empty())
        cudaMemcpy(d_cancerTypes, &cancerTypes[0], nTumors*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_gtDataMatrix, gtDataMatrix, numColT*numRowT*sizeof(bool), cudaMemcpyHostToDevice);
    cudaMemcpy(d_geDataMatrix, geDataMatrix, numColE*numRowE*sizeof(bool), cudaMemcpyHostToDevice);

    //Loop through tumors
    for(int i = rowStart; i < rowEnd; i++)
    {
        if (i % 50 == 0)
            printf("TDIC processed %d tumors.\n", i);
        
        //GET  get gt/ge indices 
        vector<int> tumorGtIndices, tumorGeIndices;
        getTumorGeneIndices(gtGeneNames[i], gtMatrixGeneNames,  tumorGtIndices);
        getTumorGeneIndices(geGeneNames[i], geMatrixGeneNames,  tumorGeIndices);
    //showing processing info on screen. 
        unsigned int nGT = tumorGtIndices.size();
        unsigned int nGE = tumorGeIndices.size();
        cout<< "Processing tumor " + tumorNames[i] + " with "<< nGT <<" Gts and "<< nGE <<" Ges \n";
        /*
        * Prepare for calling TDIC/PanCanTDIC
        */
 
        //define output rumorPosterorMatrix, this is the result of TDIC/PanCanTDIC function and will pass in to function and get results back to main
        vector<float> tumorPosteriorMatrix(nGT*nGE,0.0);

        if (!gtcFilePath.empty())
        {
            PanCanTDIC(*panCanGtMatrix, *geMatrix, globalDriverMap, tumorGtIndices, tumorGeIndices, v0, tumorPosteriorMatrix, numCanTypes, d_cancerTypes, d_gtDataMatrix, d_geDataMatrix, tumorCanTypes[i]);
        }
        else
            TDIC(*gtMatrix, *geMatrix, globalDriverMap, tumorGtIndices, tumorGeIndices, v0, tumorPosteriorMatrix, d_gtDataMatrix, d_geDataMatrix);

        //Write the contents of tumorPosteriorMatrix into a csv file
        TDIC_Output(tumorPosteriorMatrix, tumorNames[i], gtGeneNames[i], geGeneNames[i], outPath);
    }    
        // free device global memory
    if (!gtcFilePath.empty())
        cudaFree(d_cancerTypes);

    cudaFree(d_gtDataMatrix);
    cudaFree(d_geDataMatrix);

    if (!gtcFilePath.empty())
        delete panCanGtMatrix;
    else                    
        delete gtMatrix;

    delete geMatrix;

    time (&t_end);
    long seconds = difftime (t_end,t_start);
    
    int hours, minutes;
 
    minutes = seconds / 60;
    hours = minutes / 60;
    
    cout <<  " Elasped time is  " << hours << " hours " << minutes%60 << " minutes " << seconds%60 << " seconds." << "\n";

    return 0;
}
        
