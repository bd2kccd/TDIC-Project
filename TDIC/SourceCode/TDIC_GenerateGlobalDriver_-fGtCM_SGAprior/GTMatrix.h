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
//    void calcLnTumorPriors(vector<int>& gtIndx, const float v0, vector<float>& lnTumorPriorAlphas);
    void load(string gtfileName, string priorfileName);
    void calcGlobalAlphaPriors(const float v0, vector<float>& allAlphaPriors);
    float* getPriorPtr(void) const {return SGApriors;};
protected:
    float* SGApriors;

};

#endif	/* GTMATRIX_H */

