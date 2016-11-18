/*
 
 */

#ifndef TDIC_GELOOP_CUDA_H
#define TDIC_GELOOP_CUDA_H


__global__ void TDIC_GeLoop(int d_nGE, int d_nGT, int d_nTumors,  int* tumorGeIndices,  int* tumorGtIndices, int* tumorGlobDriverIndices, bool* gtDataMatrix, bool* geDataMatrix, float* lntumorMutPriors, float* tumorPosteriorMatrix);

#endif /* TDIC_GELOOP_CUDA_H */

