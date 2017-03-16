
    
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

//void PanCanTDIC(PanCanGTMatrix&, TDIMatrix&, map<string, string> & , const int, const string outPath, const float v0, int nCanType, int*& d_cancerTypes,  int*& d_gtDataMatrix, int*& d_geDataMatrix);
void PanCanTDIC(PanCanGTMatrix&, TDIMatrix&, map<string, string> & , const int, const string outPath, const float v0, int nCanType, int*& d_cancerTypes,  bool*& d_gtDataMatrix, bool*& d_geDataMatrix);
 
#endif	/* PanCanTDIC_H */

