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
    GTMatrix(string fileName);
    //GTMatrix(string fileName);
    virtual ~GTMatrix();
//    void calcLnTumorPriors(vector<int>& gtIndx, const float v0, vector<float>& lnTumorPriorAlphas);
    void load(string fileName);
    void calcGlobalAlphaPriors(const float v0, vector<float>& allAlphaPriors);
protected:


};

#endif	/* GTMATRIX_H */

