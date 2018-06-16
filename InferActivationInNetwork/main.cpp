/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: XIM33
 *
 * Created on April 9, 2018, 3:17 PM
 */
 
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
#include "Data.h"
#include <stdio.h> 
#include <math.h>
using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
    time_t t_start,t_end;
    time (&t_start);
    string phosphFile, intervFile, edgeFile, outPath;
    int hasOpt;
    while((hasOpt = getopt(argc, argv, "hp:i:e:o:")) != -1)
    {
        switch(hasOpt)
        {
            case 'p':
                phosphFile = optarg;
                break;
                
            case 'i':
                intervFile = optarg;
                break;

            case 'e':
                edgeFile = optarg;
                break;

            case 'o':
                outPath = optarg;
                break;

            case 'h':
                cerr << "Usage: InferNet -p phosphFile -i intervFile -e edgeFile -o outPath \n";
                exit(1);
                break;
                
            case '?':
                if(optopt == 'p')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;
                }
                else if(optopt == 'i')
                {
                  cout << "Option -" << optopt << " requires an argument.\n";
                  return 0;          
                }
                else if(optopt == 'e')
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
                cerr << "Usage: InferNet -p phosphFile -i intervFile -e edgeFile -o outPath \n";
                abort();
        }
    }
    
    //in construct, readin files and build network
    Data data(phosphFile, intervFile, edgeFile);

    cout << "Calculate CPT for each node.." << "\n";
    data.calCPTofEachNode();
        
    double prevJointProb = 0.0;
    int count = 0;  
    while (true){
        count ++;
    //calculate the converge parameter
        cout << "Calculate joint probability of all nodes.." << "\n";
        double jointProb = data.calJointProbOfAllNodes();
        float diff = fabs(jointProb - prevJointProb);
        cout << "Joint probability difference is " << diff << "\n";
        if (diff < 0.001 || count > 1000)  {
//        if (count > 5){
            data.outputActivatMatrix(outPath);
            data.outputCombinedMatrix(outPath);
            data.outputJointProb(outPath);
            data.outputRandomNumMatrix(outPath);
                    
            break;
        }
        else{
            prevJointProb = jointProb;
            //calculate CPT and lookup in CPT to get the inferState of each protein for each case
            cout << "Infer activation ..." << "\n";
            data.inferActivation();
            
           
        }
     }
    cout << "Number of infer iterations is " << count << "\n";
    
    time (&t_end);
    long seconds = difftime (t_end,t_start);
    
    int hours, minutes;
 
    minutes = seconds / 60;
    hours = minutes / 60;
    
    cout <<  " Elasped time is  " << hours << " hours " << minutes%60 << " minutes " << seconds%60 << " seconds." << "\n";

}



