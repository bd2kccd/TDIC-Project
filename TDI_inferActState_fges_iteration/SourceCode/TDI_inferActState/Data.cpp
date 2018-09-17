#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <string>
#include <float.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include "Data.h"



using namespace std;

Data::Data(string MatrixFileName, string SGAFileName, string DEGFileName, string edgeFileName) {

    cout << "Read in SGA matrix.. " << SGAFileName << "\n";
    readinMatrix(SGAFileName, 'S');

    if (DEGFileName != "") {
        cout << "Read in DriverActState matrix.. " << MatrixFileName << "\n";
        readinMatrix(MatrixFileName, 'M');
        cout << "Read in DEG matrix.. " << DEGFileName << "\n";
        readinMatrix(DEGFileName, 'M');
    } else {
        cout << "Read in complete matrix.. " << MatrixFileName << "\n";
        readinMatrix(MatrixFileName, 'M');
    }

    cout << "Read in edge list file " << edgeFileName << "\n";
    readinEdges(edgeFileName);

    cout << "build network ..." << "\n";
    buildNetwork();

    srand(3);

}

//Data::Data(const Data& orig) {
//}

Data::~Data() {
    for (int i = 0; i < nodeList.size(); i++) {
        delete nodeList[i];
    }
}

void Data::readinMatrix(string fileName, char type) {
    stringstream ss;
    string line;
    ifstream inFileStream;
    vector<int*> matrixAsVec;
    int nCol = 0, nRow = 0;


    FILE * pFile;
    pFile = fopen(fileName.c_str(), "r");

    if (pFile == NULL) {
        cout << "Fail to open the file " << fileName << " \n";
        exit(1);
    } else {
        fclose(pFile);
    }

    try {
        inFileStream.open(fileName.c_str());
        if ((inFileStream.rdstate() & ifstream::failbit) != 0) {
            cerr << "Error opening file when loading " + fileName + ", quit.\n";
            inFileStream.close();
            exit(EXIT_FAILURE);
        }
    } catch (...) {
        cerr << "Fail to open file " << fileName;

    }


    getline(inFileStream, line);
    if (!line.empty() && line[line.size() - 1] == '\r')
        line.erase(line.size() - 1);
    string tmp;
    stringstream firstSS(line);

    while (getline(firstSS, tmp, ',')) {
        if (type == 'M') {
            nodeNames.push_back(tmp);
        }
        nCol++;
    }

    while (getline(inFileStream, line)) {
        stringstream anotherss(line);
        string tmp;
        int curCol = 0;
        matrixAsVec.push_back(new int[nCol]());
        while (getline(anotherss, tmp, ',')) {
            matrixAsVec[nRow][curCol] = atoi(tmp.c_str());
            curCol++;
        }
        nRow++;
    }

    inFileStream.close();

    for (int c = 0; c < nCol; c++) {
        for (int r = 0; r < nRow; r++) {

            if (type == 'M') {
                combinedMatrix.push_back(matrixAsVec[r][c]);
            } else { // if (type == 'S') 
                intervMatrix.push_back(matrixAsVec[r][c]);
            }
        }
    }
    if (type == 'S') {
        nDriver = nCol;
        nTumor = nRow;
    }
    for (int i = 0; i < matrixAsVec.size(); i++)
        delete [] matrixAsVec[i];
}

void Data::readinEdges(string edgeFileName) {
    ifstream inFileStream;
    string line;
    vector<string> fields;


    FILE * pFile;
    pFile = fopen(edgeFileName.c_str(), "r");
    if (pFile == NULL) {
        cout << "Fail to open the file " << edgeFileName << " \n";
        exit(1);
    } else {
        fclose(pFile);
    }

    try {
        inFileStream.open(edgeFileName.c_str());

        while (getline(inFileStream, line)) {
            if (!line.empty() && line[line.size() - 1] == '\r')
                line.erase(line.size() - 1);
            edges.push_back(line);
        }
        inFileStream.close();
    } catch (ifstream::failure e) {
        cout << "Fail to open file " << edgeFileName;

    }
}

void Data::buildNetwork() {
    
    nNode = nodeNames.size();
    Node * node;
    for (int i = 0; i < nNode; i++) {
        string name = nodeNames[i];
        if (i < nDriver)
            node = new Node(name, 'S', i);
        else
            node = new Node(name, 'D', i);
        nodeList.push_back(node);
        nodeListMap.insert(std::pair<string, int>(name, i));
    }
////for test purpose
//    for (int i = 0; i < nDriver; i++) {
//        nodeList[i]->addToParents(i + nDriver);
//    }
//    for (int i = nDriver; i < 2 * nDriver; i++) {
//        nodeList[i]->addToChildren(i - nDriver);
//    }
//    
////test end
    for (int i = 0; i < edges.size(); i++) {

        string edge = edges[i];
        int p = edge.find(',');

        string headname = edge.substr(0, p);
        string tailname = edge.substr(p + 1);

        if (nodeListMap.find(headname) == nodeListMap.end()) {
            cout << headname << " is not found in edge " << edge << "\n";
            continue;
        }
        if (nodeListMap.find(tailname) == nodeListMap.end()) {
            cout << tailname << " is not found in edge " << edge << "\n";
            continue;
        }


        int headNodeIndex = nodeListMap[headname];
        int tailNodeIndex = nodeListMap[tailname];

        nodeList[headNodeIndex]->addToChildren(tailNodeIndex);
        nodeList[tailNodeIndex]->addToParents(headNodeIndex);

    }


    
    
}



double Data::getCPT1valueOfANode(int Aindex, int tumorID) {
    Node* Anode = nodeList[Aindex];
    vector<int>& parentIndex = Anode ->getParents();
    int numOfParents = parentIndex.size();
    vector<int> parentValues;
    for (int p = 0; p < numOfParents; p++) {
        parentValues.push_back(combinedMatrix[nTumor * parentIndex[p] + tumorID]);
    }

    vector<double>& CPT1 = Anode->getCPT();
    int lookupIndex = 0;
    for (int p = numOfParents - 1; p >= 0; p--) {
        lookupIndex += pow(2, p) * parentValues[p];
    }
    return CPT1[lookupIndex];

}

void Data::getCPTvalueOfOneChildOfANode(int childIndex, int Aindex, int tumorID, vector<double>& AchildCPT) {

    Node* childNode = nodeList[childIndex];
    int childValue = combinedMatrix [ nTumor * childIndex + tumorID];

    vector<int>& parentIndex = childNode->getParents();
    int numOfParents = parentIndex.size();
    vector<int> parentValues0, parentValues1;
    for (int p = 0; p < numOfParents; p++) {
        if (parentIndex[p] == Aindex) {
            parentValues1.push_back(1);
            parentValues0.push_back(0);
        } else {
            parentValues1.push_back(combinedMatrix[nTumor * parentIndex[p] + tumorID]);
            parentValues0.push_back(combinedMatrix[nTumor * parentIndex[p] + tumorID]);
        }
    }

    vector<double>& CPT1 = childNode->getCPT();
    int lookupIndex1 = 0;
    int lookupIndex0 = 0;
    for (int p = numOfParents - 1; p >= 0; p--) {
        lookupIndex1 += pow(2, p) * parentValues1[p];
        lookupIndex0 += pow(2, p) * parentValues0[p];
    }

    if (childValue == 1) {
        AchildCPT[1] = CPT1[lookupIndex1];
        AchildCPT[0] = CPT1[lookupIndex0];
    } else {
//        //test start
//        cout << CPT1.size() << "\n";
//        for (int k = 0 ; k < CPT1.size(); k++){
//            cout << CPT1[k] << ",";
//        }
//        cout << '\n' << lookupIndex1 << "," << lookupIndex0 << "\n";
//        //test end
        AchildCPT[1] = 1 - CPT1[lookupIndex1];
        AchildCPT[0] = 1 - CPT1[lookupIndex0];
    }

}

void Data::calCPTofEachNode() {
    for (int i = 0; i < nNode; i++) {
        nodeList[i]->calculateCPT(combinedMatrix, nTumor);
    }
}

void Data::getInferState() {
    float inferValue;

    for (int p = 0; p < nDriver; p++){
//        cout << "inferring driver " << p << " ....\n";
        for (int c = 0; c < nTumor; c++) {
            
            inferValue = inferActivation(p, c);

            float r = static_cast<float> (rand()) / static_cast<float> (RAND_MAX);

            int oriValue = combinedMatrix[p * nTumor + c];
            if (inferValue > r)
                combinedMatrix[p * nTumor + c] = 1;
            else
                combinedMatrix[p * nTumor + c] = 0;

            int newValue = combinedMatrix[p * nTumor + c];
            int change = newValue - oriValue;
            if (change != 0) {
                nodeList[p]->updateACPT(change, combinedMatrix, nTumor, c);
                vector<int>& childrenIndex = nodeList[p]->getChildren();
                for (int i = 0; i < childrenIndex.size(); i++) {
                    int idx = childrenIndex[i];
                    nodeList[idx]->updateChildCPT(change, p, combinedMatrix, nTumor, c);
                }
            }
        }
    }
}

double Data::inferActivation(int p, int c) {

    activatMatrix.clear();

    double inferValue = 0.0;

    int intervValue = intervMatrix[p * nTumor + c];

    if (intervValue == 1) {
        inferValue = 1;
    }
    else {
        inferValue = calculateProbA(p, c);
    }

    activatMatrix.push_back(inferValue);

    return inferValue;

}

double Data::calculateProbA(int p, int c) {

    double inferValue = 0.0;
    double pA1 = getCPT1valueOfANode(p, c);
    double pA0 = 1 - pA1;

    vector<double> AchildrenCPT1, AchildrenCPT0;
    vector<int>& Achildren = nodeList[p]->getChildren();

    for (int i = 0; i < Achildren.size(); i++) {
        vector<double> AchildCPT(2, 0.0);
        getCPTvalueOfOneChildOfANode(Achildren[i], p, c, AchildCPT);
        AchildrenCPT1.push_back(AchildCPT[1]);
        AchildrenCPT0.push_back(AchildCPT[0]);
    }

    double lgPA1 = 0, lgPA0 = 0;
    int TpA1 = 1, TpA0 = 1;

    if (pA1 == 0)
        TpA1 = 0;
    else
        lgPA1 = log(pA1);

    if (pA0 == 0)
        TpA0 = 0;
    else
        lgPA0 = log(pA0);

    for (int i = 0; i < Achildren.size(); i++) {

        if (AchildrenCPT1[i] == 0)
            TpA1 = 0;
        else
            lgPA1 += log(AchildrenCPT1[i]);

        if (AchildrenCPT0[i] == 0)
            TpA0 = 0;
        else
            lgPA0 += log(AchildrenCPT0[i]);
    }


    if (TpA0 == 0 && TpA1 == 0)
        inferValue = 0.5;
    else if (TpA0 == 0)
        inferValue = 1;
    else if (TpA1 == 0)
        inferValue = 0;
    else
        inferValue = 1 / (1 + exp(lgPA0 - lgPA1));

    return inferValue;

}

double Data::calJointProbOfAllNodes() {
    double logJointProb = 0.0;
    for (int p = 0; p < nNode; p++) {
        for (int c = 0; c < nTumor; c++) {
            int Avalue = combinedMatrix[p * nTumor + c];
            double pA = 0.0;
            if (Avalue == 1) {
                pA = getCPT1valueOfANode(p, c);
            } else {
                pA = 1 - getCPT1valueOfANode(p, c);
            }
            //for test purpose
            if (pA == 0.0) 
                cout << "PA == 0 in jont prob calculation\n";
            //for test end
            logJointProb += log(pA);
        }
    }
    //for test
    cout << "logJointProb \n";
    cout << logJointProb << "\n";
    jointProbs.push_back(logJointProb);
    return logJointProb;
}

void Data::outputCombinedMatrix(string outPath) {
    string outFileName;

    if (*outPath.rbegin() != '/')
        outPath = outPath + "/";

    outFileName = outPath + "combinedMatrix" + ".csv";

    FILE * pFile;
    pFile = fopen(outFileName.c_str(), "w");
    if (pFile == NULL) {
        cout << "Fail to open the file " << outFileName << " \n";
        exit(1);
    } else {
        fclose(pFile);
    }


    ofstream outFile;
    try {
        outFile.open(outFileName.c_str(), ios::out);
    } catch (ofstream::failure e) {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
    }

    for (int i = 0; i < nNode; i++) {
        outFile << nodeNames[i] ;
        if (i < nNode - 1) {
            outFile << ",";
        }
    }
    outFile << "\n";

    for (int c = 0; c < nTumor; c++) {
        for (int p = 0; p < nNode; p++) {
            outFile << combinedMatrix[p * nTumor + c];
            if (p < nNode - 1) {
                outFile << ",";
            }
        }
        outFile << "\n";
    }
    outFile.close();

}

void Data::outputActivatMatrix(string outPath) {
    string outFileName;

    if (*outPath.rbegin() != '/')
        outPath = outPath + "/";

    outFileName = outPath + "activationMatrix" + ".csv";

    FILE * pFile;
    pFile = fopen(outFileName.c_str(), "w");
    if (pFile == NULL) {
        cout << "Fail to open the file " << outFileName << " \n";
        exit(1);
    } else {
        fclose(pFile);
    }

    ofstream outFile;
    try {
        outFile.open(outFileName.c_str(), ios::out);
    } catch (ofstream::failure e) {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
    }

    for (int i = 0; i < nDriver; i++) {
        outFile << nodeNames[i];
        if (i < nDriver - 1) {
            outFile << ",";
        }
    }
    outFile << "\n";

    for (int c = 0; c < nTumor; c++) {
        for (int p = 0; p < nDriver; p++) {
            outFile << combinedMatrix[p * nTumor + c];
            if (p < nDriver - 1) {
                outFile << ",";
            }
        }
        outFile << "\n";
    }
    outFile.close();
}

void Data::outputJointProb(string outPath) {
    string outFileName;
    string outFileName2;
    if (*outPath.rbegin() != '/')
        outPath = outPath + "/";

    outFileName = outPath + "jointProbs" + ".csv";
    outFileName2 = outPath + "jointProbsAll.csv";
    
    FILE * pFile;
    FILE * pFile2;
    pFile = fopen(outFileName.c_str(), "w");
    if (pFile == NULL) {
        cout << "Fail to open the file " << outFileName << " \n";
        exit(1);
    } else {
        fclose(pFile);
    }

    //    pFile2 = fopen(outFileName2.c_str(), "w");
    //    if (pFile == NULL) {
    //        cout << "Fail to open the file " << outFileName2 << " \n";
    //        exit(1);
    //    } else {
    //        fclose(pFile);
    //    }
    

    ofstream outFile;
    try {
        outFile.open(outFileName.c_str(), ios::out);
    } catch (ofstream::failure e) {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
    }

    ofstream outFile2;
    try {
        outFile2.open(outFileName2.c_str(), ios::app);
    } catch (ofstream::failure e) {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
    }
    

    for (int i = 0; i < jointProbs.size(); i++) {
        outFile << jointProbs[i];
        if (i < jointProbs.size() - 1) {
            outFile << ",";
        }
    }
    outFile << "\n";
    outFile.close();
    
    cout << jointProbs[jointProbs.size() - 1];
    outFile2 << "," << jointProbs[jointProbs.size() - 1];
    outFile2.close();
}

void Data::outputCPT(string outPath) {
    string outFileName;
    if (*outPath.rbegin() != '/')
        outPath = outPath + "/";

    outFileName = outPath + "CPT" + ".txt";

    FILE * pFile;
    pFile = fopen(outFileName.c_str(), "w");
    if (pFile == NULL) {
        cout << "Fail to open the file " << outFileName << " \n";
        exit(1);
    } else {
        fclose(pFile);
    }


    ofstream outFile;
    try {
        outFile.open(outFileName.c_str(), ios::out);
    } catch (ofstream::failure e) {
        cout << "Exception opening output file. Please ensure you have an existing directory for file.\n";
    }




    for (int i = 0; i < nodeList.size(); i++) {
        Node* nodeA = nodeList[i];
        vector<int> parentA = nodeA->getParents();
        vector<double> &CPT1 = nodeA->getCPT();

        outFile << nodeA->getName() << "\n";

        outFile << "Parents: ";
        for (int j = parentA.size() - 1; j >= 0; j--) {
            outFile << nodeList[parentA[j]]->getType() << "_" << nodeList[parentA[j]]->getName();
            if (j > 0) {
                outFile << ",";
            }
        }
        outFile << "\n";
        outFile << "Combination: ";


        for (int j = 0; j < pow(2, parentA.size()); j++) {

            vector<int> combin_arr;
            convertToCombination(combin_arr, j);
            int combin_num = 0;
            for (int k = 0; k < combin_arr.size(); k++) {
                combin_num = combin_num * 10 + combin_arr[k];
            }
            outFile << setfill('0') << setw(parentA.size()) << combin_num;
            if (j < pow(2, parentA.size()) - 1) {
                outFile << ",";
            }
        }

        outFile << "\n";

        outFile << "CPT1: ";
        for (int j = 0; j < CPT1.size(); j++) {
            outFile << CPT1[j];
            if (j < CPT1.size() - 1) {
                outFile << ",";
            }
        }
        outFile << "\n";
        outFile << "CPT0: ";
        for (int j = 0; j < CPT1.size(); j++) {
            outFile << 1 - CPT1[j];
            if (j < CPT1.size() - 1) {
                outFile << ",";
            }
        }
        outFile << "\n\n";

    }

    outFile.close();

}

void Data::convertToCombination(vector<int> &combin_arr, unsigned int n) {

    if (n / 2 != 0) {
        convertToCombination(combin_arr, n / 2);
    }
    combin_arr.push_back(n % 2);

}