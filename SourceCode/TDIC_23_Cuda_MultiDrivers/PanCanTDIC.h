
    
using namespace std;
#include <map>
#include <vector>
#include <string>

//#include "TDIMatrix.h"
//#include "GTMatrix.h"
#include "PanCanGTMatrix.h"
#include "TDIC.h"
//#include "PanCanTDIC_GeLoop_Cuda.h"


#ifndef PanCanTDIC_H
#define	PanCanTDIC_H

void PanCanTDIC(PanCanGTMatrix&, TDIMatrix&, map<string, vector<string> >& , const int, const string outPath, const float v0, int numCanType, int*&,  bool*&, bool*&);
 
#endif	/* PanCanTDIC_H */

