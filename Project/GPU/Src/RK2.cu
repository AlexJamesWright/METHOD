#include "RK2.h"
#include <omp.h>
#include <iostream>
#include <cstdio>

#include "nvToolsExtCuda.h"
#include "nvToolsExtCudaRt.h"


// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

void RK2::step(double * cons, double * prims, double * aux, double dt)
{
  nvtxEventAttributes_t eventAttrib = {0};
  eventAttrib.version = NVTX_VERSION;
  eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib.colorType = NVTX_COLOR_ARGB;
  eventAttrib.color = 0;
  eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
  eventAttrib.message.ascii = "RK2 step";
  nvtxRangeId_t profile_id1 = nvtxRangeStartEx(&eventAttrib);

  // Syntax
  Data * d(this->data);

  // Get timestep
  if (dt <= 0) (dt=d->dt);

  // Need some work arrays
  double *p1cons, *p1prims, *p1aux, *args1, *args2;

  int Ntot(d->Nx * d->Ny * d->Nz);

  cudaHostAlloc((void **)&p1cons, sizeof(double) * Ntot * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&p1prims, sizeof(double) * Ntot * d->Nprims,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&p1aux, sizeof(double) * Ntot * d->Naux,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&args1, sizeof(double) * Ntot * d->Ncons,
                cudaHostAllocPortable);
  cudaHostAlloc((void **)&args2, sizeof(double) * Ntot * d->Ncons,
                cudaHostAllocPortable);

  nvtxEventAttributes_t eventAttrib_2 = {0};
  eventAttrib_2.version = NVTX_VERSION;
  eventAttrib_2.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib_2.colorType = NVTX_COLOR_ARGB;
  eventAttrib_2.color = 1;
  eventAttrib_2.messageType = NVTX_MESSAGE_TYPE_ASCII;
  eventAttrib_2.message.ascii = "cons2prims";
  nvtxRangeId_t profile_id2 = nvtxRangeStartEx(&eventAttrib_2);

  // Cons2prims conversion for p1 estimate stage requires old values to start
  // the rootfind
  #pragma omp parallel for
  for (int i=0; i < d->Nx; i++) {
    #pragma omp parallel for
    for (int j=0; j < d->Ny; j++) {
      #pragma omp parallel for
      for (int k=0; k < d->Nz; k++) {
        #pragma omp parallel for
        for (int var=0; var < d->Naux; var++) {
          p1aux[ID(var, i, j, k)] = aux[ID(var, i, j, k)];
        }
        #pragma omp parallel for
        for (int var=0; var < d->Nprims; var++) {
          p1prims[ID(var, i, j, k)] = prims[ID(var, i, j, k)];
        }
      }
    }
  }
  nvtxRangeEnd(profile_id2);

  nvtxRangeId_t profile_id3 = nvtxRangeStartA("fluxMethod");
  // Get first approximation of flux contribution
  this->fluxMethod->F(cons, prims, aux, d->f, args1);
  nvtxRangeEnd(profile_id3);

  nvtxEventAttributes_t eventAttrib_4 = {0};
  eventAttrib_4.version = NVTX_VERSION;
  eventAttrib_4.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib_4.colorType = NVTX_COLOR_ARGB;
  eventAttrib_4.color = 1;
  eventAttrib_4.messageType = NVTX_MESSAGE_TYPE_ASCII;
  eventAttrib_4.message.ascii = "cons2prims 2";
  nvtxRangeId_t profile_id4 = nvtxRangeStartEx(&eventAttrib_4);

  // First stage approximation
  #pragma omp parallel for
   for (int var=0; var < d->Ncons; var++) {
     #pragma omp parallel for
     for (int i=0; i < d->Nx; i++) {
       #pragma omp parallel for
       for (int j=0; j < d->Ny; j++) {
         #pragma omp parallel for
         for (int k=0; k < d->Nz; k++) {
           p1cons[ID(var, i, j, k)] = cons[ID(var, i, j, k)] - dt * args1[ID(var, i, j, k)];
         }
       }
     }
   }
   // Apply boundary conditions and get primitive and aux vars for p1
   try {
     this->model->getPrimitiveVars(p1cons, p1prims, p1aux);
   }
   catch (const std::exception& e) {
     printf("RK2 (stage 1) raises exception with following message:\n%s\n", e.what());
     throw e;
   }
   nvtxRangeEnd(profile_id4);

   nvtxRangeId_t profile_id5 = nvtxRangeStartA("bcs");
   this->bcs->apply(p1cons, p1prims, p1aux);
   nvtxRangeEnd(profile_id5);

  nvtxEventAttributes_t eventAttrib_6 = {0};
  eventAttrib_6.version = NVTX_VERSION;
  eventAttrib_6.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib_6.colorType = NVTX_COLOR_ARGB;
  eventAttrib_6.color = 1;
  eventAttrib_6.messageType = NVTX_MESSAGE_TYPE_ASCII;
  eventAttrib_6.message.ascii = "fluxMethod 2";
  nvtxRangeId_t profile_id6 = nvtxRangeStartEx(&eventAttrib_6);

   // Get second approximation of flux contribution
   this->fluxMethod->F(p1cons, p1prims, p1aux, d->f, args2);
   nvtxRangeEnd(profile_id6);

   nvtxRangeId_t profile_id7 = nvtxRangeStartA("construct solution");
   // Construct solution
   #pragma omp parallel for
   for (int var=0; var < d->Ncons; var++) {
     #pragma omp parallel for
     for (int i=0; i < d->Nx; i++) {
       #pragma omp parallel for
       for (int j=0; j < d->Ny; j++) {
         #pragma omp parallel for
         for (int k=0; k < d->Nz; k++) {
           cons[ID(var, i, j, k)] = 0.5 * (cons[ID(var, i, j, k)] + p1cons[ID(var, i, j, k)] -
                                       dt * args2[ID(var, i, j, k)]);
         }
       }
     }
   }
   nvtxRangeEnd(profile_id7);

  nvtxEventAttributes_t eventAttrib_8 = {0};
  eventAttrib_8.version = NVTX_VERSION;
  eventAttrib_8.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib_8.colorType = NVTX_COLOR_ARGB;
  eventAttrib_8.color = 1;
  eventAttrib_8.messageType = NVTX_MESSAGE_TYPE_ASCII;
  eventAttrib_8.message.ascii = "get prims";
  nvtxRangeId_t profile_id8 = nvtxRangeStartEx(&eventAttrib_8);

   // Determine new prim and aux variables
   try {
     this->model->getPrimitiveVars(cons, prims, aux);
   }
   catch (const std::exception& e) {
     printf("RK2 (corrector) raises exception with following message:\n%s\n", e.what());
     throw e;
   }
   nvtxRangeEnd(profile_id8);

   nvtxRangeId_t profile_id9 = nvtxRangeStartA("bcs 2");
   // Apply boundary conditions
   this->bcs->apply(cons, prims, aux);
   nvtxRangeEnd(profile_id9);

   // Free arrays
   cudaFreeHost(p1cons);
   cudaFreeHost(p1prims);
   cudaFreeHost(p1aux);
   cudaFreeHost(args1);
   cudaFreeHost(args2);

   nvtxRangeEnd(profile_id1);
}
