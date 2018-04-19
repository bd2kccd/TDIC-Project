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


using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {
//    time_t t_start,t_end;
//    time (&t_start);
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
    
    //in construct, readin file and build network
    Data data(phosphFile, intervFile, edgeFile);
    double prevJointProb = 0.0;
    while (true){
        //calculate the CPT for each node
        data.calCPTofEachNode();
        //calculate the converge parameter
        double jointProb = data.calJointProbOfAllNodes();
        if (jointProb < prevJointProb){
//            outputResults();
            break;
        }
        else{
            //calculate CPT and lookup in CPT to get the inferState of each protein for each case
            data.inferActivation();
            //generate threshold to cut activation table and replace to combinedMatrix
            prevJointProb = jointProb;
            break;
        }
     }

    
//    time (&t_end);
//    long seconds = difftime (t_end,t_start);
//    
//    int hours, minutes;
// 
//    minutes = seconds / 60;
//    hours = minutes / 60;
//    
//    cout <<  " Elasped time is  " << hours << " hours " << minutes%60 << " minutes " << seconds%60 << " seconds." << "\n";

}



