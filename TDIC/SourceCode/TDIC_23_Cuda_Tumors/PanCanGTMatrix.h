

#ifndef PanCanGTMATRIX_H
#define PanCanGTMATRIX_H

#include <string>
#include<vector>

//#include "TDIMatrix.h"
#include "GTMatrix.h"

class PanCanGTMatrix : public GTMatrix {
public:
    PanCanGTMatrix();
    PanCanGTMatrix(string fileName);
    virtual ~PanCanGTMatrix();
    void load(string fileName);
    
    vector<int>& getCanTypes(void) {return canTypes;};
    int getCanTypeByTumorId(int);
//    int getCanTypeByTumorName(string);
    
//private:
    vector<int> canTypes;  //For GPU version, canTypes needs to be passed into device kernel function, so changes to public
};

#endif /* PanCanGTMATRIX_H */

