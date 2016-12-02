
#ifndef PANCANTDIC_GELOOP_CUDA_H
#define PANCANTDIC_GELOOP_CUDA_H

__device__ float logSum(float lnx, float lny);
__device__ float calcFscore(float gt1,  float gt1ge1, float gt1ge0, float gt0, float gt0ge1, float gt0ge0 );
__device__ float calcA0Fscore(float gt1,  float gt1ge1, float gt1ge0, float gt0, float gt0ge1, float gt0ge0 );

 __global__ void PanCanTDIC_GeLoop(int nGE, int nGT, int nTumors, int curCanType, int* tumorGeIndices,  int* tumorGtIndices, 
         int* tumorGlobDriverIndices, int* canTypes, bool* gtDataMatrix, bool* geDataMatrix, float* lntumorMutPriors, 
         float* tumorPosteriorMatrix, float* tumorPosteriorMatrixC );

#endif /* PANCANTDIC_GELOOP_CUDA_H */

