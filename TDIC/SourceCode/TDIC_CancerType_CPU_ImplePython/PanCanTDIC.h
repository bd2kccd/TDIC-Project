
    
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

//#ifdef	__cplusplus
//extern "C" {
//#endif

   
    //Function declarations
    void PanCanTDIC(PanCanGTMatrix&, TDIMatrix&, map<string, string> & , const int, const string outPath, const double v0, int numCanType);
 
    
    double calcPanCanFscore(const bool* gtDataMatrix, const bool* geDataMatrix, const int nTumors, const double curGtPrior,
        const int curGtRowStart, const int DriverRowStart, const int curGeRowStart, const int curGtIndx, vector<int> canTypes);

    
    


//#ifdef	__cplusplus
//}
//#endif

#endif	/* PanCanTDIC_H */

