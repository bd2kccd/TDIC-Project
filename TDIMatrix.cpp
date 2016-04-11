/* 
 * File:   TDIMatrix.cpp
 * Author: Kevin Lu
 * 
 * Created on April 25, 2015, 6:13 PM
 */
using namespace std;

#include "TDIMatrix.h"
#include <algorithm>    // std::find
#include <vector>
#include <fstream> 
#include <iostream>
#include <string>
#include <string.h>
#include <sstream>

//using namespace std;

TDIMatrix::TDIMatrix(){
    mat = NULL;
    // other initialization stuff
}

TDIMatrix::TDIMatrix(string fileName){
    this->load(fileName);    
    
}

TDIMatrix::TDIMatrix(const TDIMatrix& orig) {
    // copy constructor not implemented.
}


TDIMatrix::~TDIMatrix() {
    if (mat) delete [] mat; 
}

string TDIMatrix::getTumorNameById(const int indx)   
{
    if (tumorNames.empty())
    {
        cerr << "Bug: attempt to query an unpopulated TDIMatrix\n"; 
        return string();
    }
    return tumorNames.at(indx);
}


/**
 * 
 * @param geneIndx index of gene of interest
 * @param tumorIndx index of tumor of interest
 * @return mutation or expresssion value of the given gene for the given tumor
 */
int TDIMatrix::valueAt(int geneIndx, int tumorIndx){
    return mat[geneIndx * nTumors + tumorIndx];
}



/**
 * 
 * @param inGeneNames list of tumors
 * @param outGeneIndices corresponding indeces for the tumors in "inGeneNames"
 */
void TDIMatrix::getTumorIndicesByNames(vector<string>& inGeneNames, vector<int>& outGeneIndices)
{  

    for(int i = 0; i < inGeneNames.size(); i++)
    {
        int count = 0;
        for(int j = 0; j < tumorNames.size(); j++)
        {
            count ++;
            if(inGeneNames[i].compare(tumorNames[j]) == 0)
            //if(s.compare(tumorNames[j]) == 0)
            {
                outGeneIndices.push_back(j);
                break;
            }
        }

        if(count == geneNames.size() && inGeneNames[i].compare(geneNames[count-1]) != 0)
        {
            cout << "Bug: calling getTumorIndicesByNames fail to look up index for  gene " << inGeneNames[i] << " does not exist\n";
        }

    }
    if (inGeneNames.size() > outGeneIndices.size())
    {
        cout << "Bug: in getTumorIndicesByNames, some gene index missed\n";
    }
}


/**
 * 
 * @param inGeneNames list of genes
 * @param outGeneIndices corresponding indeces for the genes in "inGeneNames"
 */
void TDIMatrix::getGeneIndicesByNames(vector<string>& inGeneNames, vector<int>& outGeneIndx){
    for(int i = 0; i < inGeneNames.size(); i++)
    {
        int count = 0;
        for(int j = 0; j < geneNames.size(); j++)
        {
            count ++;
            if(inGeneNames[i].compare(geneNames[j]) == 0)
            {
                outGeneIndx.push_back(j);
                break;
            }
        }
    }
    if(inGeneNames.size() != outGeneIndx.size())
    {
        cout << "Bug: In getGeneIndicesByNames, try to look up gene name that does not exist\n";
        cout << "Input size: " << inGeneNames.size() << " " << "output size: " << outGeneIndx.size() << "\n";
    }
}

/**
 * 
 * @param tumorIndx list of tumor indeces
 * @param outTumorName list of tumors corresponding to the list of indices in "tumorIndx"
 */
void TDIMatrix::getTumorNamesByIndices(vector<int>& tumorIndx, vector<string>& outTumorName) {
    for(int i = 0; i < tumorIndx.size(); i++)
    {
        outTumorName.push_back(tumorNames[tumorIndx[i]]);
    }
}

/**
 * 
 * @param geneIndx list of gene indeces
 * @param outGeneName list of genes corresponding to the list of indices in "geneIndx"
 */
void TDIMatrix::getGeneNamesByIndices(vector<int>& geneIndx, vector<string>& outGeneName){
    //cout << "In getGeneNamesByIIndices\n";
    //cout << "length of gene names: " << geneNames.size() << "\n";
    for(int i = 0; i < geneIndx.size(); i++)
    {
        //cout << "geneIndex for i= " << i << ": " << geneIndx[i] << "\n";
        // cout << "geneName: " << geneNames[geneIndx[i]] << "\n";
        outGeneName.push_back(geneNames[geneIndx[i]]);
    }    
}

/**
 * 
 * @param fileName filepath of CSV file to be converted into a TDIMatrix object
 */
void TDIMatrix::load(string fileName)
{
    std::stringstream ss;
    std::string line;
    ifstream inFileStream;   
    vector<int*> matrixAsVec;
    int nCol = 0, nRow = 0;

    try
    {
        inFileStream.open(fileName.c_str());
        if ( (inFileStream.rdstate() & std::ifstream::failbit ) != 0 )
        {
            std::cerr << "Error opening file when loading TDI matrix, quit.\n";
            inFileStream.close();
            exit(EXIT_FAILURE);
        }
    }
    catch (...) 
    { //std::ifstream::failure e
        cerr << "Fail to open file " << fileName;
        
    }

    //cout << "Opened GE file ok.\n";
            
    // read in first line to figure out the number of columns and column Names;
    getline(inFileStream, line);
    line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

    string tmp;
    stringstream firstSS(line);

    bool firstColFlag = true;
    while (getline(firstSS, tmp, ',' ))
    {
        if(firstColFlag)
        {
            firstColFlag = false;
            continue;
        }
        geneNames.push_back(tmp);
        //geneIndxMap.insert(std::pair<string, int>(tmp, nCol));
        nCol++;
    }

    while (getline(inFileStream, line))
    {
        firstColFlag = true; 

        stringstream anotherss(line);
        string tmp;
        int curCol = 0;
        while (getline(anotherss, tmp, ','))
        {
            if(firstColFlag)
            {
                firstColFlag = false;
                tumorNames.push_back(tmp);
                matrixAsVec.push_back(new int[nCol]);
                continue;
            }

            matrixAsVec[nRow][curCol] = atoi(tmp.c_str());
            curCol++;
        }        
        nRow++;		
    }

    inFileStream.close();          
    
    // transform the vector of inter arrays to a consecutive array so that 
    // computation can be done efficiently
    mat = new int[nCol * nRow];
    int indx = 0;    
    for (int g = 0; g < nCol; g++){
        for (int t = 0; t < nRow; t++){
            mat[indx++] = matrixAsVec[t][g];
        }
    }
    
    nTumors = tumorNames.size();
    nGenes = geneNames.size();
    
    // free temporary matrix
    for (int i = 0; i < matrixAsVec.size(); i++) 
        delete [] matrixAsVec[i];
}


/**
 * 
 * @param geneID index of gene
 * @param tumorIndices vector of indices representing tumors that have a mutation in the given gene
 */
void TDIMatrix::findTumorsWithOnesPerGene(int geneID, vector<int>& tumorIndices)
{
    for (int t = 0; t < nTumors; t++){
        if (mat[nTumors * geneID + t] == 1) tumorIndices.push_back(t);
    }    
}

/**
 * 
 * @param tumorID index of tumor
 * @param geneIndices vector of indices representing genes that are mutated in the given tumor
 */
 //should this return A0 if it is a 1?
void TDIMatrix::findGeneWithOnesInTumor(int tumorID, vector<int>& geneIndices){
    for (int g = 0; g < nGenes; g ++ ){
        if (mat[g * nTumors + tumorID] == 1) geneIndices.push_back(g);
    }
}

bool TDIMatrix::writeToCSV(string outFilePath)
{
    ofstream outFile;
    try
    {
        outFile.open(outFilePath.c_str());
    }
    catch(ofstream::failure e)
    {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
        return false;
    }
    
    //start writing CSV representation of TDIMatrix
    
    //write column headers
    for(int i = 0; i < geneNames.size(); i++)
    {
        outFile << "," << geneNames[i];
    }
    outFile << "\n";
    
    for(int i = 0; i < tumorNames.size(); i++)
    {
        outFile << tumorNames[i];
        for(int j = 0; j < geneNames.size(); j++)
        {
            outFile << "," << mat[j * nTumors + i];
        }
        outFile << "\n";
    }
    
    outFile.close();
    return true;
}
