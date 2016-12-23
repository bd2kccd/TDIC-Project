#include <cstdlib>
#include <unistd.h>
#include <getopt.h>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <ctype.h>
#include <fstream>
#include <omp.h>
#include "TDIC.h"
#include "PanCanTDIC.h"

//#include "TDIMatrix.h"
//#include "GTMatrix.h"
//#include "PanCanGTMatrix.h"


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
                cerr << "Usage: TDIC -i inputDataFile -c inputPancanGtMatrix [-f inputGtMatrix] -d inputGeMatrix -g inputGlobDriverDictionary -o pathForOutputResults -i inputTumorFile -t canTypeCoding\n";
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
    
    cout << "Reading global driver file.\n";
    map<string, string> globalDriverMap;
    parseGlobDriverDict(globalDriverPath, globalDriverMap);

    
    /*
     * read in input data file and get tumor gt ge gene indices
     */
    
    //prepare for calling TDIC_Load function
     
    //first get TDIC gt/ge gene names
    vector<string> gtMatrixGeneNames, geMatrixGeneNames;
    if (!gtcFilePath.empty())
        gtMatrixGeneNames = panCanGtMatrix->getGeneNames();
    else
        gtMatrixGeneNames = gtMatrix->getGeneNames();
    geMatrixGeneNames = geMatrix->getGeneNames();
    
    //then get input data gt/ge names
    vector<vector<string> > gtGeneNames, geGeneNames;
    vector<string> tumorNames;
    vector<int> tumorCanTypes;
    TDIC_Load(inputDataFile, gtGeneNames, geGeneNames, tumorNames, tumorCanTypes, cancerTypeTable);
    
    int nTumors =  tumorNames.size();
    
    if(rowStart < 0)
        rowStart = 0;
    if(rowEnd < 1 || rowEnd > nTumors)
        rowEnd = nTumors;

	if(rowStart > rowEnd)
    {
        cout << "Given rowEnd index is smaller than given rowStart. Exiting out.\n";
        exit(1);
    }
        
    
//    cout<< "gtGeneName.size(), geGeneNames.size()"<<gtGeneNames.size()<<","<<geGeneNames.size()<<"\n";
    for (int i = rowStart; i< rowEnd; i++){
        // get gt/ge indices
        vector<int> gtGeneIndices, geGeneIndices;
        getTumorGeneIndices(gtGeneNames[i], gtMatrixGeneNames,  gtGeneIndices);
        getTumorGeneIndices(geGeneNames[i], geMatrixGeneNames,  geGeneIndices);
   
        /*
         * Prepare for calling TDIC/PanCanTDIC
         */

        //showing processing info on screen. 
        int nGT = gtGeneIndices.size();
        int nGE = geGeneIndices.size();
        cout << "Processing tumor " << tumorNames[i] << " with " << nGT << " GAs, and " << nGE << " GEs" << "\n";  

        //define output rumorPosterorMatrix, this is the result of TDIC/PanCanTDIC function and will pass in to function and get results back to main
        vector<float> tumorPosteriorMatrix(nGT*nGE,0.0);

        if (!gtcFilePath.empty())
        {
            PanCanTDIC(*panCanGtMatrix, *geMatrix, globalDriverMap, gtGeneIndices, geGeneIndices, v0, tumorPosteriorMatrix, tumorCanTypes[i]);
        }
        else
            TDIC(*gtMatrix, *geMatrix, globalDriverMap, gtGeneIndices, geGeneIndices, v0, tumorPosteriorMatrix);

        //Write the contents of tumorPosteriorMatrix into a csv file
        TDIC_Output(tumorPosteriorMatrix, tumorNames[i], gtGeneNames[i], geGeneNames[i], outPath);
    }
    if (!gtcFilePath.empty())
        delete panCanGtMatrix;
    else                    
        delete gtMatrix;

    delete geMatrix;
 
    return 0;
}

