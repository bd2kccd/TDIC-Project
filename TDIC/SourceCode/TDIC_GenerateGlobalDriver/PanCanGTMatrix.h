

#ifndef PanCanGTMATRIX_H
#define PanCanGTMATRIX_H

#include <string>
#include<vector>

#include "TDIMatrix.h"
#include "GTMatrix.h"

class PanCanGTMatrix : public GTMatrix {
public:
    PanCanGTMatrix();
    PanCanGTMatrix(string fileName);
    virtual ~PanCanGTMatrix();
    void load(string fileName);
    
    vector<int>& getCanTypes(void) {return canTypes;};
    int getCanTypeByTumorId(int);
    int getCanTypeByTumorName(string);
    int getNumCanType();
private:
    vector<int> canTypes;
};

#endif /* PanCanGTMATRIX_H */

