
    
using namespace std;
#include <map>
#include <vector>
#include <string>

//#include "TDIMatrix.h"
//#include "GTMatrix.h"
#include "PanCanGTMatrix.h"
#include "TDIC.h"
//#include "PanCanTDIC_GeLoop.h"


#ifndef PanCanTDIC_H
#define	PanCanTDIC_H

//void PanCanTDIC(PanCanGTMatrix&, TDIMatrix&, map<string, string> & , const int, const string outPath, const float v0, int*&,  int*&, int*&);
 void PanCanTDIC(PanCanGTMatrix& gtMatrix, TDIMatrix& geMatrix, map<string, 
        string>& mapGlobDrivers, vector<int>  tumorGtIndices, vector<int> tumorGeIndices, const float v0, vector<float>& tumorPosteriorMatrix,  
         int numCanTypes, int  *&d_cancerTypes,  bool *&d_gtDataMatrix, bool *&d_geDataMatrix, int tumorCanType);

#endif	/* PanCanTDIC_H */

