
    
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

#ifdef	__cplusplus
extern "C" {
#endif

   
    //Function declarations
    void PanCanTDIC(PanCanGTMatrix&, TDIMatrix&, map<string, string> & , const int, const string outPath, const float v0);
 
    
    
    
    


#ifdef	__cplusplus
}
#endif

#endif	/* PanCanTDIC_H */

