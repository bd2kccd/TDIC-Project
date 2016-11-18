
    
using namespace std;
#include <map>
#include <vector>
#include <string>

//#include "TDIMatrix.h"
//#include "GTMatrix.h"
#include "PanCanGTMatrix.h"
#include "TDIC.h"


#ifndef PanCanTDIC_H
#define	PanCanTDIC_H

    //Function declarations
    void PanCanTDIC(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        vector<string> >& mapGlobDrivers, vector<int>  tumorGtIndices, vector<int> tumorGeIndices, const float v0, vector<float>& tumorPosteriorMatrix);
 
#endif	/* PanCanTDIC_H */

