/*
 * This class Inherited from class GTMatrix. It handles the PanCan GTMatrix 
 * with Cancer Type at the last column in the input data file
 */

/* 
 * File:   PanCanGTMatrix.cpp
 * Author: XIM33
 * 
 * Created on March 16, 2016, 12:51 PM
 */
#include <math.h>

//#include "TDIMatrix.h"
//#include "GTMatrix.h"

using namespace std;

#include <algorithm>    // std::find
#include <vector>
#include <fstream> 
#include <iostream>
#include <string>
#include <string.h>
#include <sstream>

#include "PanCanGTMatrix.h"


PanCanGTMatrix::PanCanGTMatrix() {

   
}
PanCanGTMatrix::PanCanGTMatrix(string fileName) 
{
    this->load(fileName);
    calcGlobalAlphaPriors();
}


PanCanGTMatrix::~PanCanGTMatrix(){
   
}





/**
 * Load GT with cancel type data from raw data (CSV file).
 * The last columns are cancel types which are loaded into canTypes vector
 * @param fileName Filepath for a CSV file representing a PanCanGT matrix.
 */
void PanCanGTMatrix::load(string fileName){
    std::stringstream ss;
    std::string line;
    ifstream inFileStream;   
    vector<bool*> matrixAsVec;
    int nCol = 0, nRow = 0;

    try{
        inFileStream.open(fileName.c_str());
        if ( (inFileStream.rdstate() & std::ifstream::failbit ) != 0 )
        {
            std::cerr << "Error opening file when loading GTMC matrix, quit.\n File Name: " << fileName;
            inFileStream.close();
            exit(EXIT_FAILURE);
        }
    }
    catch (...) { //std::ifstream::failure e
        cerr << "Fail to open file " << fileName;
        
    } 
            
    // read in first line to figure out the number of columns and column Names;
    getline(inFileStream, line);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    
    string tmp;
    stringstream firstSS (line);

    bool firstColFlag = true;
    geneNames.push_back("A0"); 
    while (getline(firstSS, tmp, ',' )){
        if(firstColFlag)
        {
                firstColFlag = false;
                continue;
        }
        geneNames.push_back(tmp);
       
        //geneIndxMap.insert(std::pair<string, int>(tmp, nCol));
        nCol ++;
    }


    
    //Check the last column to see if the column name something like 'Cancel Type'
    // Cancel type column name can be 'Cancer Type' or 'Cancel Types' or 'Can Type' or 'Can Types' . 
    //There can be no or more space between each word. Each character in the word is not case sensitive, can be lower or upper case.
    
    
    
    string colName;
    colName = geneNames.back();
    //Convert to lower case
    transform(colName.begin(), colName.end(), colName.begin(), ::tolower);
    //remove all spaces
    colName.erase(remove(colName.begin(), colName.end(), ' '), colName.end());

    if (colName == "cancertype" || colName == "cancertypes" || colName == "cantype" ||colName == "cantypes"
            || colName == "cantypecode" || colName == "cancertypecode" || colName == "cantypescode" || colName == "cancertypescode")//correct column name
    {
        //last column is cancel type not gene name, so popback, nCol minus 1
        geneNames.pop_back();
        nCol--;
    }
    else//incorrect column name. 
    {
        cout << "Input GTMatrix after -c parameter or PanCan GTMatrix last column name is not correct. It must be something like \"Cancer Type\" ";
        exit(1);
    }
    while (getline(inFileStream, line)){
        firstColFlag = true; 

        stringstream anotherss (line);
        string tmp;
        int curCol = 0;
        int i_canType = 0;
        while (getline(anotherss, tmp, ',')){
                if(firstColFlag)
                {
                    firstColFlag = false;
                    tumorNames.push_back(tmp);
                    matrixAsVec.push_back(new bool[nCol]());
                    continue;
                }
                if(curCol != nCol)
                {
                    matrixAsVec[nRow][curCol] = atoi(tmp.c_str());
                    curCol++;
                }
                else
                {
                    i_canType = atoi(tmp.c_str());
                    if(i_canType <=0)
                    {
                       cout << "Error: Cancer type value is incorrect. It must be an integer bigger than zero.\n" ;
                       exit(1);
                    }
                            
                    canTypes.push_back(i_canType) ;
                }
        }
        nRow++;     
    }
    inFileStream.close();            
    

    
    // transform the vector of inter arrays to a consecutive array so that 
    // computation can be done efficiently
    int indx = 0;
    mat = new bool[(nCol+1) * nRow ]();
    for (int i = 0; i < nRow; i++)
    {
        mat[indx++] = 1;
    }
    
    for (int g = 0; g < nCol; g++){
        for (int t = 0; t < nRow; t++){
            mat[indx++] = matrixAsVec[t][g];
        }
    }
    
    nTumors = tumorNames.size();
    nGenes = geneNames.size();
    
    // free temporary matrix
    for (int i = 0; i < matrixAsVec.size(); i++) 
        delete[] matrixAsVec[i];
}

//int PanCanGTMatrix::getCanTypeByTumorId(int tumorId){
//    if (!canTypes.empty())
//        return(canTypes.at(tumorId));
//    else
//    {
//        cerr << "Cancer type vector is empty\n";
//        return -1;
//    }
//}
//int PanCanGTMatrix::getCanTypeByTumorName(string tumorName){
//    for (int i = 0; i < tumorNames.size(); i++)
//    {
//        if (tumorNames.at(i).compare(tumorName) == 0)
//            return canTypes.at(i);
//    }
//    //if rumorName not found, then pop out the error message
//    cerr << "Tumor name is not found in getCanTypeByTumorName function";
//    return -1;
//}
