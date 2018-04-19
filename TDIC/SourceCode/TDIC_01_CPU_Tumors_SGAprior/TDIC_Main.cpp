// This version read in tumor-SGA prior matrix. Prior is not calculated from GtMatrix anymore.
// When loading GtMatrix, at the same time load the prior matrix as well

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
    
    double v0 = 0.1;
    int hasOpt;
    string gtFilePath, gtcFilePath, globalDriverPath, degFilePath, outPath, inputDataFile, cancerTypeTable;
    int rowStart = -1;
    int rowEnd = -1;

    time_t t_start,t_end;
    time (&t_start);

    while((hasOpt = getopt(argc, argv, "hi:f:d:g:o:s:e:c:t:?")) != -1)
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
                    exit(0);
                }
                break;
            
            case 'e':
                rowEnd = atoi(optarg);
                if(rowEnd < 0)
                {
                    cout << "rowEnd given is less than zero. Exiting out.\n";
                    exit(0);                    
                }
                break;

            case 'h':
                cerr << "Usage: TDIC -i inputData -p tumorPriorFile -f inputGtMatrix -d inputGeMatrix -g inputGlobDriverDictionary -o pathForOutputResults \n";
                exit(0);
                break;
                
            case '?':
                if(optopt == 'p')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;
                }
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
                else if(optopt == 't')
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
                cerr << "Usage: TDIC -i inputDataFile -p tumorPriorFile -c inputPancanGtMatrix [-f inputGtMatrix] -d inputGeMatrix -g inputGlobDriverDictionary -o pathForOutputResults -i inputTumorFile -t canTypeCoding\n";
                abort();
        }
    }
     
    
    if (inputDataFile.empty() )
    {
        cerr << "Must have input data file\n";
        exit(0);
    }
    else if (!gtFilePath.empty() && !gtcFilePath.empty() )//input both GTMatrix and PanCanGTMatrix
    {
        //gtFileePath and gtcFilePath can not both exist, either process GTMatrix or PanCanGTMatrix
        cerr << "Can not input both GtMatrix and PanCanGtMatrix\n";
        exit(0);            
    }
    else if (gtFilePath.empty() && gtcFilePath.empty() )//both gtFileePath and gtcFilePath not exist
    {
        
        cerr << "Must input GtMatrix or PanCanGtMatrix \n";
        exit(0);            
    }    
    else if (globalDriverPath.empty())
    {
        cerr << "Must input global driver file\n";
        exit(1);
    }
    else if (cancerTypeTable.empty())
    {
        cerr << "Must input cancer type coding file\n";
        exit(1);
    }
    
    
    GTMatrix* gtMatrix;
    PanCanGTMatrix* panCanGtMatrix;
    
    if (!gtFilePath.empty()) //input GTMatrix
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
    vector<vector<string> > gtNamesAll, geNamesAll;
    vector<string> tumorNames;
    vector<int> tumorCanTypes;
    vector<vector<double> > gtPriorsAll;
    TDIC_Load(inputDataFile, gtNamesAll, gtPriorsAll, geNamesAll, tumorNames, tumorCanTypes, cancerTypeTable);
    
    int nTumors =  tumorNames.size();
    
    if(rowStart < 0)
        rowStart = 0;
    if(rowEnd < 1 || rowEnd > nTumors)
        rowEnd = nTumors;

	if(rowStart >= rowEnd)
    {
        cout << "Given rowEnd index must bigger than given rowStart. Exiting out.\n";
        exit(1);
    }
        
    for (int i = rowStart; i< rowEnd; i++){
        // get gt/ge indices
        vector<int> gtIndicesTumor, geIndicesTumor;
        vector<double> gtPriorsTumor = gtPriorsAll[i];
        getTumorGeneIndices(gtNamesAll[i], gtMatrixGeneNames,  gtIndicesTumor);
        getTumorGeneIndices(geNamesAll[i], geMatrixGeneNames,  geIndicesTumor);
   
        /*
         * Prepare for calling TDIC/PanCanTDIC
         */

        //showing processing info on screen. 
        int nGT = gtIndicesTumor.size();
        int nGE = geIndicesTumor.size();
        cout << "Processing tumor " << tumorNames[i] << " with " << nGT << " GAs, and " << nGE << " GEs" << "\n";  

        //define output rumorPosterorMatrix, this is the result of TDIC/PanCanTDIC function and will pass in to function and get results back to main
        vector<double> tumorPosteriorMatrix(nGT*nGE,0.0);

        if (!gtcFilePath.empty())
        {
            PanCanTDIC(*panCanGtMatrix, *geMatrix, globalDriverMap, gtIndicesTumor, geIndicesTumor, v0, tumorPosteriorMatrix, tumorCanTypes[i]);
        }
        else
            TDIC(*gtMatrix, *geMatrix, globalDriverMap, gtIndicesTumor,gtPriorsTumor, geIndicesTumor, v0, tumorPosteriorMatrix);

        //Write the contents of tumorPosteriorMatrix into a csv file
        TDIC_Output(tumorPosteriorMatrix, tumorNames[i], gtNamesAll[i], geNamesAll[i], outPath);
    }
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

