// In this version, global prior is not calculated on the fly. It is generated previously
//and loaded into TDI by -p parameter

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
#include "GeneGlobDriver.h"
#include "GenePanCanGlobDriver.h"
//#include "TDIMatrix.h"
//#include "GTMatrix.h"
//#include "PanCanGTMatrix.h"


using namespace std;

int main(int argc, char** argv) {
    
    time_t t_start,t_end;
    time (&t_start);    
    
    
    int hasOpt;
    int nTumors;
    float v0 = 0.01;//this is the default value if not input from command line. 
    
    
    GTMatrix* gtMatrix;
    PanCanGTMatrix* panCanGtMatrix;
//    string gtFilePath, gtcFilePath, globalDriverPath, degFilePath, outPath;
    string gtFilePath, gtcFilePath, degFilePath, outPath,strv0, priorFilePath;

    while((hasOpt = getopt(argc, argv, "h:f:d:v:o:c:p:?")) != -1)
    {
        switch(hasOpt)
        {
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

            case 'v':
                strv0 = optarg ;
                v0 = atof(strv0.c_str());
                break;
            
            case 'o':
                outPath = optarg;
                break;

            case 'h':
                cerr << "Usage: TDIC -p globalPrior -f inputGaMatrix -d inputGeMatrix -o pathForOutputResults [ -v0 constant v0]\n";
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
                cerr << "Usage: TDIC -p globalPrior -f inputGaMatrix -d inputGeMatrix -v constant v0 -o pathForOutputResults\n";
                abort();
        }
    }
    
    cout<<"v0="<< v0 <<"\n";

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
        cout << "Reading GT matrix: " << gtFilePath << " and \n  global prior file: " << priorFilePath << "\n";
        gtMatrix = new GTMatrix(gtFilePath,priorFilePath);
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
   


    if (!gtFilePath.empty())//process GTMatrix
    {

        GeneGlobDriver(*gtMatrix, *geMatrix, outPath, v0);
        delete gtMatrix;
    }
    else//process PanCanGTMatrix
    {
  
        GenePanCanGlobDriver(*panCanGtMatrix, *geMatrix, outPath, v0);
        delete panCanGtMatrix;
    }
      
    delete geMatrix;  
    
    time (&t_end);
    long seconds = difftime (t_end,t_start);
    
    int hours, minutes;
 
    minutes = seconds / 60;
    hours = minutes / 60;
    
    cout <<  " Elasped time is  " << hours << " hours " << minutes%60 << " minutes " << seconds%60 << " seconds.";

    

    return 0;
}

