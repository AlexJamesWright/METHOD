#include "saveData.h"
#include <cstdlib>
#include <cstdio>

using namespace std;

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

void SaveData::saveAll(bool timeSeries)
{
  // Determine the directory to write files to
  if (!timeSeries) {
    strcpy(dir, "Data/Final");
    app[0]=0;
  }
  else {
    strcpy(dir, "Data/TimeSeries");
    sprintf(app, "%d", Nouts++);
  }


  this->saveCons();
  this->savePrims();
  this->saveAux();
  this->saveConsts();
}

void SaveData::saveCons()
{
  FILE * f;

  char fname[60];
  strcpy(fname, dir);
  strcat(fname, "/Conserved/cons");
  strcat(fname, app);
  strcat(fname, ".dat");

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


  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          fprintf(f, "%.16f ", d->cons[ID(var, i, j, k)]);
        }
        fprintf(f, "\n");
      }
    }
  }
  fclose(f);

}


void SaveData::savePrims()
{
  FILE * f;
  char fname[60];
  strcpy(fname, dir);
  strcat(fname, "/Primitive/prims");
  strcat(fname, app);
  strcat(fname, ".dat");  f = fopen(fname, "w");

  // Ensure file is open
  if (f == NULL) {
    printf("Error: could not open 'prims.dat' for writing.\n");
    exit(1);
  }

  // File is open, write data
  fprintf(f, "prims = ");
  for (int i(0); i < d->Nprims-1; i++) fprintf(f, "%s, ", d->primsLabels[i].c_str());
  fprintf(f, "%s\n", d->primsLabels[d->Nprims-1].c_str());
  for (int var(0); var < d->Nprims; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          fprintf(f, "%.16f ", d->prims[ID(var, i, j, k)]);
        }
        fprintf(f, "\n");
      }
    }
  }

  fclose(f);

}

void SaveData::saveAux()
{
  FILE * f;
  char fname[60];
  strcpy(fname, dir);
  strcat(fname, "/Auxilliary/aux");
  strcat(fname, app);
  strcat(fname, ".dat");  f = fopen(fname, "w");

  // Ensure file is open
  if (f == NULL) {
    printf("Error: could not open 'aux.dat' for writing.\n");
    exit(1);
  }

  // File is open, write data
  fprintf(f, "aux = ");
  for (int i(0); i < d->Naux-1; i++) fprintf(f, "%s, ", d->auxLabels[i].c_str());
  fprintf(f, "%s\n", d->auxLabels[d->Naux-1].c_str());
  for (int var(0); var < d->Naux; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          fprintf(f, "%.16f ", d->aux[ID(var, i, j, k)]);
        }
        fprintf(f, "\n");
      }
    }
  }

  fclose(f);

}


void SaveData::saveConsts()
{
  FILE * f;
  char fname[60];
  strcpy(fname, dir);
  strcat(fname, "/Constants/constants");
  strcat(fname, app);
  strcat(fname, ".dat");  f = fopen(fname, "w");

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
