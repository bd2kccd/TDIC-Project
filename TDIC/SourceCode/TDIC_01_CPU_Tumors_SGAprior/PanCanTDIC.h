
    
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
        string>& mapGlobDrivers, vector<int>  tumorGtIndices, vector<int> tumorGeIndices, const double v0, vector<double>& tumorPosteriorMatrix, int tumorCanType);
 
#endif	/* PanCanTDIC_H */

