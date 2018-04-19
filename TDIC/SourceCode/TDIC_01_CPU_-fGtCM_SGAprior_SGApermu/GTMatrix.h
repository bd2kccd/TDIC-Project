/* 
 * File:   GTMatrix.h
 * Author: kevin
 *
 * Created on May 21, 2015, 9:36 AM
 */

#ifndef GTMATRIX_H
#define	GTMATRIX_H
#include "TDIMatrix.h"

class GTMatrix : public TDIMatrix{
public:
    GTMatrix();
    GTMatrix(string gtfileName, string priorfileName);
    //GTMatrix(string fileName);
    virtual ~GTMatrix();
    void calcLnTumorPriors(vector<int>& gtIndx, const float v0, vector<float>& lnTumorPriorAlphas);
    void load(string gtfileName, string priorfileName);
    float* getPriorPtr(void) const {return mat_prior;};
    
 protected:
    void calcGlobalAlphaPriors();
    float* allAlphaPriors;

    float* mat_prior; //  mat_prior has the same data structure as the mat, the difference is that it stores prior instead of gt value
                // we represent the matrix as a consecutive liner array 
                // with  gene-by-tumor.  This way, it is more efficient to 
                // collect statistics from tumors because we loop through tumors

};

#endif	/* GTMATRIX_H */

