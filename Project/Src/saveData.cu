#include "saveData.h"
#include <cstdlib>
#include <cstdio>
#include <string>

void SaveData::saveAll()
{
  this->saveCons();
  this->savePrims();
  this->saveAux();
  this->saveConsts();
}

void SaveData::saveCons()
{
  FILE * f;
  f = fopen("Data/conserved.dat", "w");

  // Ensure file is open
  if (f == NULL) {
    printf("Error: could not open 'conserved.dat' for writing.\n");
    exit(1);
  }

  // File is open, write data
  fprintf(f, "cons = ");
  for (int i(0); i < d->Ncons-1; i++) fprintf(f, "%s, ", d->consLabels[i].c_str());
  fprintf(f, "%s\n", d->consLabels[d->Ncons-1].c_str());
  for (int var(0); var < d->Ncons; var++) {
    for (int i(0); i < d->Nx; i++) {
      for (int j(0); j < d->Ny; j++) {
        for (int k(0); k < d->Nz; k++) {
          fprintf(f, "%f ", d->cons[d->id(var, i, j, k)]);
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
  f = fopen("Data/primitive.dat", "w");

  // Ensure file is open
  if (f == NULL) {
    printf("Error: could not open 'primitive.dat' for writing.\n");
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
          fprintf(f, "%f ", d->prims[d->id(var, i, j, k)]);
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
  f = fopen("Data/auxilliary.dat", "w");

  // Ensure file is open
  if (f == NULL) {
    printf("Error: could not open 'auxilliary.dat' for writing.\n");
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
          fprintf(f, "%f ", d->aux[d->id(var, i, j, k)]);
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
  f = fopen("Data/constants.dat", "w");

  // Ensure file is open
  if (f == NULL) {
    printf("Error: could not open 'constants.dat' for writing.\n");
    exit(1);
  }

  fprintf(f, "constants = nx, ny, nz, Nx, Ny, Nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime, cfl, Ng, gamma, sigma, ");
  fprintf(f, "Ncons, Nprims, Naux, cp, dt, t, dx, dy, dz\n");
  fprintf(f, "%d %d %d %d %d %d %f %f %f %f %f %f %f %f %d %f %f %d %d %d %f %f %f %f %f %f\n",
          d->nx, d->ny, d->nz, d->Nx, d->Ny, d->Nz, d->xmin, d->xmax, d->ymin, d->ymax, d->zmin, d->zmax, d->endTime, d->cfl, d->Ng,
          d->gamma, d->sigma, d->Ncons, d->Nprims, d->Naux, d->cp, d->dt, d->t, d->dx, d->dy, d->dz);

  fclose(f);

}
