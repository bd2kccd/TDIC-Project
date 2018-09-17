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
#include <stdio.h> 
#include <math.h>


#include "Data.h"

using namespace std;

int main(int argc, char** argv) {
    time_t t_start,t_end;
    time (&t_start);
    string MatrixFile="", DEGFile=""; //matrixFile can be completeMatrix or DriverActState
    string SGAFile, edgeFile, outPath;
    int iterNum = 2000;

    int hasOpt;
    while((hasOpt = getopt(argc, argv, "hm:s:e:o:d:x:")) != -1)
    {
        switch(hasOpt)
        {
            case 'm':
                MatrixFile = optarg;
                break;
                
            case 's':
                SGAFile = optarg;
                break;

            case 'e':
                edgeFile = optarg;
                break;

            case 'd':
                DEGFile = optarg;
                break;
             
            case 'x':
                iterNum = atoi(optarg);
                break;
                
            case 'o':
                outPath = optarg;
                break;

            case 'h':
                cerr << "Usage: InferNet -m completeFile -s SGAFile -e edgeFile -o outPath or\n";
                cerr << "Usage: InferNet -m drivrActStateFile -d DEGfile -s SGAFile -e edge File -o outPath \n";
                exit(1);
                break;
                
            case '?':
                if(optopt == 'm')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;
                }
                else if(optopt == 's')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;          
                }
                else if(optopt == 'e')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;          
                }

                else if(optopt == 'd')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;          
                }
                else if(optopt == 'x')
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
                cerr << "Usage: InferNet -m completeFile -s SGAFile -e edgeFile -o outPath or\n";
                cerr << "Usage: InferNet -m drivrActStateFile -d DEGfile -s SGAFile -e edge File -o outPath \n";
                abort();
        }
    }
    
    
    Data data(MatrixFile, SGAFile, DEGFile,edgeFile);

    cout << "Calculate CPT for each node.." << "\n";
    data.calCPTofEachNode();
        
    double prevJointProb = 0.0;
    int count = 0;  
    while (true){
        count ++;
    
        cout << "Calculate joint probability of all nodes.." << "\n";
        double jointProb = data.calJointProbOfAllNodes();
        float diff = fabs(jointProb - prevJointProb);
        cout << "Joint probability difference is " << diff << "\n";
        if (diff < 0.001 || count >= iterNum)  {
            data.outputActivatMatrix(outPath);
            data.outputCombinedMatrix(outPath);
            data.outputJointProb(outPath);
            data.outputCPT(outPath);
            
                    
            break;
        }
        else{
            prevJointProb = jointProb;
            cout << "Inferring activation count " << count << "\n";
            data.getInferState();
            
           
        }
     }
    cout << "Number of total iterations is " << count << "\n";
    
    time (&t_end);
    long seconds = difftime (t_end,t_start);
    
    int hours, minutes;
 
    minutes = seconds / 60;
    hours = minutes / 60;
    
    cout <<  " Elasped time is  " << hours << " hours " << minutes%60 << " minutes " << seconds%60 << " seconds." << "\n";

}



