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
#include <time.h>
#include "TDIC.h"
#include "PanCanTDIC.h"
//#include "TDIMatrix.h"
//#include "GTMatrix.h"
//#include "PanCanGTMatrix.h"


using namespace std;

int main(int argc, char** argv) {
    
    // parse arguments 
    //extern char *optarg;
    //extern int optind, opterr, optopt;
    
    
    time_t t_start,t_end;
    time (&t_start);
    
    int rowStart = -1;
    int rowEnd = -1;
    int hasOpt;
    int nTumors;
    GTMatrix* gtMatrix;
    PanCanGTMatrix* panCanGtMatrix;
    string gtFilePath, gtcFilePath, globalDriverPath, degFilePath, outPath, priorFilePath,gtGene;

    while((hasOpt = getopt(argc, argv, "hs:e:f:d:g:o:c:p:t:")) != -1)
    {
        switch(hasOpt)
        {
            case 't':
                gtGene = optarg;
                break;
                
            case 'p':
                priorFilePath = optarg;
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

            case 'o':
                outPath = optarg;
                break;

            case 'h':
                cerr << "Usage: TDIC -f inputGaMatrix -d inputGeMatrix -g inputGlobDriverDictionary -o pathForOutputResults [-s rowStart index -e rowEnd index]\n";
                exit(1);
                break;
                
            case '?':
                if(optopt == 'f')
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
                cerr << "Usage: TDIC -f inputGaMatrix -d inputGeMatrix -g inputGlobDriverDictionary -o pathForOutputResults\n";
                abort();
        }
    }

    if (priorFilePath.empty())
    {
        cerr <<"Must input tumor-SGA prior matrix";
        exit(1);
    }
    
    if (!gtFilePath.empty() && !gtcFilePath.empty() )//input both GTMatrix and PanCanGTMatrix
    {
        //gtFileePath and gtcFilePath can not both exist, either process GTMatrix or PanCanGTMatrix
        cerr << "Can not input both GtMatrix and PanCanGtMatrix\n";
        exit(1);            
    }
    else if (gtFilePath.empty() && gtcFilePath.empty() )
     {
        //both gtFileePath and gtcFilePath not exist
        cerr << "Must input GtMatrix or PanCanGtMatrix \n";
        exit(1);            
    }    
    else if (!gtFilePath.empty()) //input GTMatrix
    {           //read in GT matrices
        cout << "Reading GT matrix: " << gtFilePath << "\n";
        gtMatrix = new GTMatrix(gtFilePath, priorFilePath);
        nTumors = gtMatrix->getNTumors();
    }
    else //input PanCanGTMatrix
    {   cout << "Reading PanCanGT matrix: " << gtcFilePath << "\n";
        panCanGtMatrix = new PanCanGTMatrix(gtcFilePath);
        nTumors = panCanGtMatrix->getNTumors();
    }   
     
 
    //read in GE matrices
       
    cout << "Reading GE matrix. " << degFilePath << "\n";
    TDIMatrix* geMatrix = new TDIMatrix(degFilePath);
   
    cout << "Reading global driver file.\n";
    map<string, string> globalDriverMap;
    parseGlobDriverDict(globalDriverPath, globalDriverMap);
    
    if(rowStart == -1)
        rowStart = 0;
    if(rowEnd == -1 || rowEnd > nTumors)
        rowEnd = nTumors;

	if(rowStart > rowEnd)
    {
        cout << "Given rowEnd index is smaller than given rowStart. Exiting out.\n";
        exit(1);
    }

    //tumorNames = gtMatrix->getTumorNames();
//    vector<int> outGlobDriverIndx;

    // check use gpu flag, if yes branch out to TCIGPU(GTmatrix, DEG, Glog)
    
    
    float v0 = 0.1;
    if (!gtFilePath.empty())//process GTMatrix
    {
        
        bool* gtDataMatrix = gtMatrix->getMatPtr();
        vector<string> gtGeneNames =  gtMatrix->getGeneNames();
        int gtPos = std::find(gtGeneNames.begin(), gtGeneNames.end(), gtGene) - gtGeneNames.begin(); 
        if(gtPos >= gtGeneNames.size()) {
            cout << "Gene " << gtGene << " not found \n";
            exit(1);
        }
//        #pragma omp parallel for
        for(int i = rowStart; i < rowEnd; i++)
        {
            // if (i % 50 == 0)
                // printf("TDIC processed %d tumors.\n", i);

            
            if (gtDataMatrix[nTumors*gtPos+i] == 1)
            {
                printf("TDIC is processing %dth tumor.\n", i);
                TDIC(*gtMatrix, *geMatrix, globalDriverMap, i, outPath, v0);
            }

	}	
        delete gtMatrix;
			        
        
//        for(int i = rowStart; i < rowEnd; i++)
//        {
//            if (i % 50 == 0)
//                printf("TDIC processed %d tumors.\n", i);
//            TDIC(*gtMatrix, *geMatrix, globalDriverMap, i, outPath, v0);
//                    
//        }
//        delete gtMatrix;
    }
    else//process PanCanGTMatrix
    {
        //get number of caner types 
        int NcanType = panCanGtMatrix->getNumCanType();
		cout << "There are " << NcanType << " cancer types. \n";
//        #pragma omp parallel for
        for(int i = rowStart; i < rowEnd; i++)
        {
            if (i % 50 == 0)
                printf("TDIC processed %d tumors.\n", i);
            PanCanTDIC(*panCanGtMatrix, *geMatrix, globalDriverMap, i, outPath, v0);
        }
        delete panCanGtMatrix;
    }

    delete geMatrix;  

    time (&t_end);
    long seconds = difftime (t_end,t_start);
    
    int hours, minutes;
 
    minutes = seconds / 60;
    hours = minutes / 60;
    
    cout <<  " Elasped time is  " << hours << " hours " << minutes%60 << " minutes " << seconds%60 << " seconds." << "\n";

    return 0;
}

