#include "parallelSaveData.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <mpi.h>

using namespace std;

// Id in a state vector that does not include ghost cells
// TODO -- Should probably just define a variable on Data that is (Nz-2*Ng or 1 if nz=0) to avoid having a copy for each dimension
#define ID_PHYS_3D(variable, idx, jdx, kdx) ((variable)*(d->Nx-(d->Ng*2))*(d->Ny-(d->Ng*2))*(d->Nz-(d->Ng*2)) + (idx)*(d->Ny-(d->Ng*2))*(d->Nz-(d->Ng*2)) + (jdx)*(d->Nz-(d->Ng*2)) + (kdx))
#define ID_PHYS_2D(variable, idx, jdx) ((variable)*(d->Nx-(d->Ng*2))*(d->Ny-(d->Ng*2)) + (idx)*(d->Ny-(d->Ng*2)) + (jdx))
#define ID_PHYS_1D(variable, idx) ((variable)*(d->Nx-(d->Ng*2)) + (idx))

#define ID_FULL_3D(variable, idx, jdx, kdx) ((variable)*(d->nx)*(d->ny)*(d->nz) + (idx)*(d->ny)*(d->nz) + (jdx)*(d->nz) + (kdx))
#define ID_FULL_2D(variable, idx, jdx) ((variable)*(d->nx)*(d->ny) + (idx)*(d->ny) + (jdx))
#define ID_FULL_1D(variable, idx) ((variable)*(d->nx) + (idx))
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

void ParallelSaveData::saveAll(bool timeSeries)
{
  // Clean directory variable
  dir[0] = '\0';
  // Determine the directory to write files to
  if (test)
    strcpy(dir, "../../");
  if (!timeSeries && strcmp(dir, "Data/Final")!=0) {
    strcat(dir, "Data/Final");
    app[0]=0;
  }
  else {
    strcat(dir, "Data/TimeSeries");
    sprintf(app, "%d", Nouts++);
  }

  // Cons
  this->saveCons();

  // Prims
  this->savePrims();

  // Aux
  this->saveAux();

  // TODO -- could gather this to proc0 like for the other state vectors but not sure if it is required
  //this->saveDomain();

  // TODO -- Nx, Ny are per process -- may need to print out a global version as well (nx, ny don't include ghost cells)
  this->saveConsts();

}

void ParallelSaveData::packStateVectorBuffer(double *buffer, double *stateVector, int nVars){
  // Prepare send buffer, which doesn't include ghost cells, by copying from local state vectors
  if (d->dims==3){
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->Nx-(d->Ng*2); i++) {
        for (int j(0); j < d->Ny-(d->Ng*2); j++) {
          for (int k(0); k < d->Nz-(d->Ng*2); k++) {
            buffer[ID_PHYS_3D(var, i, j, k)] = stateVector[ID(var, i + d->Ng, j + d->Ng, k + d->Ng)];
          }
        }
      }
    }
  } else if (d->dims==2){
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->Nx-(d->Ng*2); i++) {
        for (int j(0); j < d->Ny-(d->Ng*2); j++) {
          buffer[ID_PHYS_2D(var, i, j)] = stateVector[ID(var, i + d->Ng, j + d->Ng, 0)];
        }
      }
    }
  } else {
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->Nx-(d->Ng*2); i++) {
        buffer[ID_PHYS_1D(var, i)] = stateVector[ID(var, i + d->Ng, 0, 0)];
      }
    }
  }
}

void ParallelSaveData::copyMasterStateVectorToFullStateVector(double *fullStateVector, double *stateVector, int nVars){
  // This requires proc0 to have xRankId=yRankId=zRankId=0
  if (d->dims==3){
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->Nx-(d->Ng*2); i++) {
        for (int j(0); j < d->Ny-(d->Ng*2); j++) {
          for (int k(0); k < d->Nz-(d->Ng*2); k++) {
            fullStateVector[ID_FULL_3D(var, i, j, k)] = stateVector[ID(var, i + d->Ng, j + d->Ng, k + d->Ng)];
          }
        }
      }
    }
  } else if (d->dims==2){
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->Nx-(d->Ng*2); i++) {
        for (int j(0); j < d->Ny-(d->Ng*2); j++) {
          //printf("nx: %d, ny: %d\n", d->nx, d->ny);
          //printf("var: %d i: %d j: %d, id: %d, id_full: %d\n", var, i, j, ID(var, i+d->Ng, j+d->Ng, 0),
                  //ID_FULL_2D(var, i, j));
          fullStateVector[ID_FULL_2D(var, i, j)] = stateVector[ID(var, i + d->Ng, j + d->Ng, 0)];
        }
      }
    }
  } else {
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->Nx-(d->Ng*2); i++) {
        fullStateVector[ID_FULL_1D(var, i)] = stateVector[ID(var, i + d->Ng, 0, 0)];
      }
    }
  }
}

void ParallelSaveData::sendStateVectorBufferToMaster(double *buffer, int numCellsSent, int rank){
   // MPI message vars
   int tag = 101;
   MPI_Status status;
   if (env->rank == rank){
       MPI_Send(buffer, numCellsSent, MPI_DOUBLE, 0, tag, env->mpiCartesianComm);
   } else if (env->rank == 0){
       MPI_Recv(buffer, numCellsSent, MPI_DOUBLE, rank, tag, env->mpiCartesianComm, &status);
   }
}

void ParallelSaveData::unpackStateVectorBuffer(double *buffer, double *stateVector, int nVars, int rank){
  // Unpack send buffer, which don't include ghost cells, into the global state vector

  // Get (x,y,z) coords of rank that sent data to proc0
  int rankCoords[3];
  int ndims = 3; // rank grid is always 3D
  MPI_Cart_coords(env->mpiCartesianComm, rank, ndims, rankCoords);

  int iOffset, jOffset, kOffset;
  iOffset = rankCoords[0] * (d->Nx - (d->Ng*2));
  if (d->dims > 1) {
      jOffset = rankCoords[1] * (d->Ny - (d->Ng*2));
  } else jOffset = 0;

  if (d->dims > 2) {
      kOffset = rankCoords[2] * (d->Nz - (d->Ng*2));
  } else kOffset = 0;

  if (d->dims==3){
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->Nx-(d->Ng*2); i++) {
        for (int j(0); j < d->Ny-(d->Ng*2); j++) {
          for (int k(0); k < d->Nz-(d->Ng*2); k++) {
            stateVector[ID_FULL_3D(var, i + iOffset, j + jOffset, k + kOffset)] = buffer[ID_PHYS_3D(var, i, j, k)];
          }
        }
      }
    }
  } else if (d->dims==2){
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->Nx-(d->Ng*2); i++) {
        for (int j(0); j < d->Ny-(d->Ng*2); j++) {
          stateVector[ID_FULL_2D(var, i + iOffset, j + jOffset)] = buffer[ID_PHYS_2D(var, i, j)];
        }
      }
    }
  } else {
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->Nx-(d->Ng*2); i++) {
        stateVector[ID_FULL_1D(var, i + iOffset)] = buffer[ID_PHYS_1D(var, i)];
      }
    }
  }
}

void ParallelSaveData::writeStateVectorToFile(FILE *f, double *fullStateVector, int nVars){
  if (d->dims==3){
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->nx; i++) {
        for (int j(0); j < d->ny; j++) {
          for (int k(0); k < d->nz; k++) {
            fprintf(f, "%.16f ", fullStateVector[ID_FULL_3D(var, i, j, k)]);
          }
          fprintf(f, "\n");
        }
      }
    }
  } else if (d->dims==2){
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->nx; i++) {
        for (int j(0); j < d->ny; j++) {
          fprintf(f, "%.16f ", fullStateVector[ID_FULL_2D(var, i, j)]);
          fprintf(f, "\n");
        }
      }
    }
  } else {
    for (int var(0); var < nVars; var++) {
      for (int i(0); i < d->nx; i++) {
        fprintf(f, "%.16f ", fullStateVector[ID_FULL_1D(var, i)]);
        fprintf(f, "\n");
      }
    }
  }
}

void ParallelSaveData::saveCons()
{
  FILE * f;

  char fname[120];
  strcpy(fname, dir);
  strcat(fname, "/Conserved/cons");
  strcat(fname, app);
  strcat(fname, ".dat\0");

  // Allocate buffers for gathering distributed state vectors onto master process
  // We do this here rather than in saveAll to allow saveCons to be called independently
  // We don't want to do this in the ParallelSaveData constructor as we don't want to use up this large
  // amount of memory until it's needed.
  int numCellsInBuffer = d->Ncons * (d->Nx-(2*d->Ng));
  if (d->dims > 1) numCellsInBuffer *= (d->Ny - (2*d->Ng));
  if (d->dims > 2) numCellsInBuffer *= (d->Nz - (2*d->Ng));
  double *buffer = (double*) malloc(numCellsInBuffer * sizeof(double));
  int numCellsInFullStateVector = numCellsInBuffer * env->nProc;
  double *fullStateVector = (double*) malloc(numCellsInFullStateVector * sizeof(double));

  // For all procs other than proc0, copy local statevector to a buffer that does not include ghost cells
  // for sending to proc0. Proc0 can copy directly from its local statevector to the fullstatevector
  if (env->rank != 0) packStateVectorBuffer(buffer, d->cons, d->Ncons);
  else copyMasterStateVectorToFullStateVector(fullStateVector, d->cons, d->Ncons);

  for (int r(1); r < env->nProc; r++){
      int numCellsSent = d->Ncons * (d->Nx-(2*d->Ng));
      if (d->dims > 1) numCellsSent *= (d->Ny-(2*d->Ng));
      if (d->dims > 2) numCellsSent *= (d->Nz-(2*d->Ng));
      sendStateVectorBufferToMaster(buffer, numCellsSent, r);
      if (env->rank == 0) unpackStateVectorBuffer(buffer, fullStateVector, d->Ncons, r);
  }

  if (env->rank == 0){
    f = fopen(fname, "w");
    // Ensure file is open
    if (f == NULL) {
      printf("Error: could not open 'cons.dat' for writing.\n");
      exit(1);
    }

    // File is open, write data
    fprintf(f, "cons = ");
    for (int i(0); i < d->Ncons-1; i++) {
      fprintf(f, "%s, ", d->consLabels[i].c_str());
    }
    fprintf(f, "%s\n", d->consLabels[d->Ncons-1].c_str());

    writeStateVectorToFile(f, fullStateVector, d->Ncons);

    fclose(f);
  }

  free(buffer);
  free(fullStateVector);
}

void ParallelSaveData::savePrims()
{
  FILE * f;
  char fname[120];
  strcpy(fname, dir);
  strcat(fname, "/Primitive/prims");
  strcat(fname, app);
  strcat(fname, ".dat\0");  f = fopen(fname, "w");

  // Allocate buffers for gathering distributed state vectors onto master process
  // We do this here rather than in saveAll to allow savePrims to be called independently
  // We don't want to do this in the ParallelSaveData constructor as we don't want to use up this large
  // amount of memory until it's needed.
  int numCellsInBuffer = d->Nprims * (d->Nx-(2*d->Ng));
  if (d->dims > 1) numCellsInBuffer *= (d->Ny - (2*d->Ng));
  if (d->dims > 2) numCellsInBuffer *= (d->Nz - (2*d->Ng));
  double *buffer = (double*) malloc(numCellsInBuffer * sizeof(double));
  int numCellsInFullStateVector = numCellsInBuffer * env->nProc;
  double *fullStateVector = (double*) malloc(numCellsInFullStateVector * sizeof(double));

  if (env->rank != 0) packStateVectorBuffer(buffer, d->prims, d->Nprims);
  else copyMasterStateVectorToFullStateVector(fullStateVector, d->prims, d->Nprims);
  for (int r(1); r < env->nProc; r++){
      int numCellsSent = d->Nprims * (d->Nx-(2*d->Ng));
      if (d->dims > 1) numCellsSent *= (d->Ny-(2*d->Ng));
      if (d->dims > 2) numCellsSent *= (d->Nz-(2*d->Ng));
      sendStateVectorBufferToMaster(buffer, numCellsSent, r);
      if (env->rank == 0) unpackStateVectorBuffer(buffer, fullStateVector, d->Nprims, r);
  }

  if (env->rank == 0){
    // Ensure file is open
    if (f == NULL) {
      printf("Error: could not open 'prims.dat' for writing.\n");
      exit(1);
    }

    // File is open, write data
    fprintf(f, "prims = ");
    for (int i(0); i < d->Nprims-1; i++) fprintf(f, "%s, ", d->primsLabels[i].c_str());
    fprintf(f, "%s\n", d->primsLabels[d->Nprims-1].c_str());

    writeStateVectorToFile(f, fullStateVector, d->Nprims);
    fclose(f);
  }

  free(buffer);
  free(fullStateVector);
}

void ParallelSaveData::saveAux()
{
  FILE * f;
  char fname[120];
  strcpy(fname, dir);
  strcat(fname, "/Auxiliary/aux");
  strcat(fname, app);
  strcat(fname, ".dat\0");  f = fopen(fname, "w");

  // Allocate buffers for gathering distributed state vectors onto master process
  // We do this here rather than in saveAll to allow saveAux to be called independently
  // We don't want to do this in the ParallelSaveData constructor as we don't want to use up this large
  // amount of memory until it's needed.
  int numCellsInBuffer = d->Naux * (d->Nx-(2*d->Ng));
  if (d->dims > 1) numCellsInBuffer *= (d->Ny - (2*d->Ng));
  if (d->dims > 2) numCellsInBuffer *= (d->Nz - (2*d->Ng));
  double *buffer = (double*) malloc(numCellsInBuffer * sizeof(double));
  int numCellsInFullStateVector = numCellsInBuffer * env->nProc;
  double *fullStateVector = (double*) malloc(numCellsInFullStateVector * sizeof(double));

  if (env->rank != 0) packStateVectorBuffer(buffer, d->aux, d->Naux);
  else copyMasterStateVectorToFullStateVector(fullStateVector, d->aux, d->Naux);
  for (int r(1); r < env->nProc; r++){
      int numCellsSent = d->Naux * (d->Nx-(2*d->Ng));
      if (d->dims > 1) numCellsSent *= (d->Ny-(2*d->Ng));
      if (d->dims > 2) numCellsSent *= (d->Nz-(2*d->Ng));
      sendStateVectorBufferToMaster(buffer, numCellsSent, r);
      if (env->rank == 0) unpackStateVectorBuffer(buffer, fullStateVector, d->Naux, r);
  }

  if (env->rank == 0){
    // Ensure file is open
    if (f == NULL) {
      printf("Error: could not open 'aux.dat' for writing.\n");
      exit(1);
    }

    // File is open, write data
    fprintf(f, "aux = ");
    for (int i(0); i < d->Naux-1; i++) fprintf(f, "%s, ", d->auxLabels[i].c_str());
    fprintf(f, "%s\n", d->auxLabels[d->Naux-1].c_str());

    writeStateVectorToFile(f, fullStateVector, d->Naux);
    fclose(f);
  }

  free(buffer);
  free(fullStateVector);

}


void ParallelSaveData::saveDomain()
{
  FILE * f;
  char fname[120];
  strcpy(fname, dir);
  strcat(fname, "/Domain/domain");
  strcat(fname, app);
  strcat(fname, ".dat\0");  f = fopen(fname, "w");

  // Ensure file is open
  if (f == NULL) {
    printf("Error: could not open 'domain.dat' for writing.\n");
    exit(1);
  }

  // File is open, write data
  for (int i(0); i < d->Nx; i++)
    fprintf(f, "%.16f ", d->x[i]);
  fprintf(f, "\n");
  for (int j(0); j < d->Ny; j++)
    fprintf(f, "%.16f ", d->y[j]);
  fprintf(f, "\n");
  for (int k(0); k < d->Nz; k++)
    fprintf(f, "%.16f ", d->z[k]);
  fprintf(f, "\n");


  fclose(f);

}


void ParallelSaveData::saveConsts()
{
  FILE * f;
  char fname[120];
  strcpy(fname, dir);
  strcat(fname, "/Constants/constants");
  strcat(fname, app);
  strcat(fname, ".dat\0");  f = fopen(fname, "w");

  if (env->rank == 0){
      // Ensure file is open
      if (f == NULL) {
        printf("Error: could not open 'constants.dat' for writing.\n");
        exit(1);
      }

      fprintf(f, "constants = nx, ny, nz, Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, cfl, Ng, gamma, sigma, ");
      fprintf(f, "Ncons, Nprims, Naux, cp, dt, t, dx, dy, dz\n");
      fprintf(f, "%d %d %d %d %d %d %.16f %.16f %.16f %.16f %.16f %.16f %.16f %.16f %d %.16f %.16f %d %d %d %.16f %.16f %.16f %.16f %.16f %.16f\n",
              d->nx, d->ny, d->nz, d->Nx, d->Ny, d->Nz, d->xmin, d->xmax, d->ymin, d->ymax, d->zmin, d->zmax, d->endTime, d->cfl, d->Ng,
              d->gamma, d->sigma, d->Ncons, d->Nprims, d->Naux, d->cp, d->dt, d->t, d->dx, d->dy, d->dz);

      fclose(f);
  }
}


void ParallelSaveData::saveVar(string variable, int num)
{
  int cpa(0); // cons=1,prims=2,aux=3
  int Nvar(0); // Variable number
  FILE * f;
  char fname[120];
  double * sendVec; // Pointer to the array to send to master and save

  // Determine which variable the user wants saved
  for (int var(0); var < d->Ncons; var++) {
    if (strcmp(d->consLabels[var].c_str(), variable.c_str()) == 0) {
      cpa=1; Nvar=var;
      break;
    }
  }

  if (!cpa) {
    for (int var(0); var < d->Nprims; var++) {
      if (strcmp(d->primsLabels[var].c_str(), variable.c_str()) == 0) {
        cpa=2; Nvar=var;
        break;
      }
    }
  }

  if (!cpa) {
    for (int var(0); var < d->Naux; var++) {
      if (strcmp(d->auxLabels[var].c_str(), variable.c_str()) == 0) {
        cpa=3; Nvar=var;
        break;
      }
    }
  }

  if (!cpa) {
    printf("Error: Could not find user specified variable '%s'\n", variable.c_str());
    exit(1);
  }

  if (cpa==1) sendVec = &d->cons[ID(Nvar, 0, 0, 0)];
  else if (cpa==2) sendVec = &d->prims[ID(Nvar, 0, 0, 0)];
  else sendVec = &d->aux[ID(Nvar, 0, 0, 0)];

  // Allocate buffers for gathering distributed state vectors onto master process
  // We do this here rather than in saveAll to allow savePrims to be called independently
  // We don't want to do this in the ParallelSaveData constructor as we don't want to use up this large
  // amount of memory until it's needed.
  int numCellsInBuffer = (d->Nx-(2*d->Ng));
  if (d->dims > 1) numCellsInBuffer *= (d->Ny - (2*d->Ng));
  if (d->dims > 2) numCellsInBuffer *= (d->Nz - (2*d->Ng));
  double *buffer = (double*) malloc(numCellsInBuffer * sizeof(double));
  int numCellsInFullStateVector = numCellsInBuffer * env->nProc;
  double *fullStateVector = (double*) malloc(numCellsInFullStateVector * sizeof(double));

  if (env->rank != 0) packStateVectorBuffer(buffer, sendVec, 1);
  else copyMasterStateVectorToFullStateVector(fullStateVector, sendVec, 1);
  for (int r(1); r < env->nProc; r++){
      int numCellsSent = 1 * (d->Nx-(2*d->Ng));
      if (d->dims > 1) numCellsSent *= (d->Ny-(2*d->Ng));
      if (d->dims > 2) numCellsSent *= (d->Nz-(2*d->Ng));
      sendStateVectorBufferToMaster(buffer, numCellsSent, r);
      if (env->rank == 0) unpackStateVectorBuffer(buffer, fullStateVector, 1, r);
  }






  if (env->rank == 0){

    // Directory
    if (this->test)
      strcpy(fname, "../../Data/TimeSeries/UserDef/");
    else
      strcpy(fname, "Data/TimeSeries/UserDef/");
    sprintf(app, "%d", Nouts);

    // Location of output file
    strcat(fname, variable.c_str());
    strcat(fname, app);
    strcat(fname, ".dat\0");
    f = fopen(fname, "w");

    // Ensure file is open
    if (f == NULL) {
      printf("Error: could not open user-defined file for writing.\n");
      exit(1);
    }

    // File is open, write data
    fprintf(f, "var = %s, t = %18.16f\n", variable.c_str(), d->t);

    writeStateVectorToFile(f, fullStateVector, 1);


    fclose(f);


    // For first output add the variables we are saving
    if (Nouts==0) {
      if (Ncount==0) {
        ofstream info;
        if (this->test)
          strcpy(fname, "../../Data/TimeSeries/UserDef/");
        else
          strcpy(fname, "Data/TimeSeries/UserDef/");
        strcat(fname, "info");
        info.open(fname);
        info << variable << endl;
        info.close();
      }
      else {
        ofstream info;
        info.open("Data/TimeSeries/UserDef/info", ios::app);
        info << variable << endl;
        info.close();
      }
    }
    Ncount++;
    // Increment if this is the last variable to save in this timestep
    if (Ncount == num) {
      Ncount = 0;
      Nouts++;
    }
  }

  free(buffer);
  free(fullStateVector);

}
