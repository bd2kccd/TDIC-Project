

#ifndef PanCanGTMATRIX_H
#define PanCanGTMATRIX_H

#include <string>
#include<vector>
#include <algorithm>

#include "TDIMatrix.h"
#include "GTMatrix.h"


class PanCanGTMatrix : public GTMatrix {
public:
    PanCanGTMatrix();
    PanCanGTMatrix(string fileName);
    virtual ~PanCanGTMatrix();
    void load(string fileName);
    
//    vector<int>& getCanTypes(void) {return canTypes;};
    int getCanTypeByTumorId(int);
//    int getCanTypeByTumorName(string);
    //get number of cancer types. (using max value of the cancer type as the number of cancer types)
    int getNumCanType(void){return *max_element(canTypes.begin(), canTypes.end());};
    
private:
    vector<int> canTypes;
};

#endif /* PanCanGTMATRIX_H */

