#ifndef CUDAERRORCHECK_H
#define CUDAERRORCHECK_H

#include <stdio.h>
#include <stdlib.h>

// Error checking for CUDA calls
inline void gpuAssert(cudaError_t code, const char *file, int line)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"\n\nGPUassert: '%s' in %s, line %d\n\n\n", cudaGetErrorString(code), file, line);
      exit(code);
   }
}

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }


#endif
