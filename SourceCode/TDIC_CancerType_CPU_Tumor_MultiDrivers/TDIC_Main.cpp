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
    string gtFilePath, gtcFilePath, globalDriverPath, degFilePath, outPath, inputDataFile ;

    while((hasOpt = getopt(argc, argv, "hi:f:d:g:o:c:?")) != -1)
    {
        switch(hasOpt)
        {
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
                cerr << "Usage: TDIC -i inputDataFile -c inputPancanGtMatrix [-f inputGtMatrix] -d inputGeMatrix -g inputGlobDriverDictionary -o pathForOutputResults\n";
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
    map<string, vector<string> > globalDriverMap;
    parseGlobDriverDict(globalDriverPath, globalDriverMap);
    
//    //**********for test****************
//    for(map<string, vector<string> >::iterator it =globalDriverMap.begin(); it != globalDriverMap.end(); ++it) {
//        cout << it->first << "\n";
//        vector<string> drivers = it->second;
//        for (int i=0; i<drivers.size(); i++)
//            cout <<drivers[i]<<",";
//        cout << "\n";
//    }
//    
//    //**********for text end***************

    
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
    vector<string> gtGeneNames, geGeneNames;
    string curTumorName;
    TDIC_Load(inputDataFile, gtGeneNames, geGeneNames, curTumorName);
    
    //last to get gt/ge indices 
    vector<int> gtGeneIndices, geGeneIndices;
    getTumorGeneIndices(gtGeneNames, gtMatrixGeneNames,  gtGeneIndices);
    getTumorGeneIndices(geGeneNames, geMatrixGeneNames,  geGeneIndices);
   
    /*
     * Prepare for calling TDIC/PanCanTDIC
     */
    
    //showing processing info on screen. 
    int nGT = gtGeneIndices.size();
    int nGE = geGeneIndices.size();
    cout << "Processing tumor " << curTumorName << " with " << nGT << " GAs, and " << nGE << " GEs" << "\n";  

    //define output rumorPosterorMatrix, this is the result of TDIC/PanCanTDIC function and will pass in to function and get results back to main
    vector<float> tumorPosteriorMatrix(nGT*nGE,0.0);
   
    if (!gtcFilePath.empty())
    {
        PanCanTDIC(*panCanGtMatrix, *geMatrix, globalDriverMap, gtGeneIndices, geGeneIndices, v0, tumorPosteriorMatrix);
    }
    else
        TDIC(*gtMatrix, *geMatrix, globalDriverMap, gtGeneIndices, geGeneIndices, v0, tumorPosteriorMatrix);

    //Write the contents of tumorPosteriorMatrix into a csv file
    TDIC_Output(tumorPosteriorMatrix, curTumorName, gtGeneNames, geGeneNames, outPath);

    if (!gtcFilePath.empty())
        delete panCanGtMatrix;
    else                    
        delete gtMatrix;

    delete geMatrix;
 
    return 0;
}

