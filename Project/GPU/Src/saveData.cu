#include "saveData.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>

using namespace std;

// Macro for getting array index
#define ID(variable, idx, jdx, kdx) ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

void SaveData::saveAll(bool timeSeries)
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

  this->saveCons();
  this->savePrims();
  this->saveAux();
  this->saveDomain();
  this->saveConsts();
}

void SaveData::saveCons()
{
  FILE * f;

  char fname[60];
  strcpy(fname, dir);
  strcat(fname, "/Conserved/cons");
  strcat(fname, app);
  strcat(fname, ".dat\0");

  f = fopen(fname, "w");
  // Ensure file is open
  printf("Writing into %s\n", fname);
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
  strcat(fname, ".dat\0");  f = fopen(fname, "w");

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
  strcat(fname, "/Auxiliary/aux");
  strcat(fname, app);
  strcat(fname, ".dat\0");  f = fopen(fname, "w");

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


void SaveData::saveDomain()
{
  FILE * f;
  char fname[60];
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


void SaveData::saveConsts()
{
  FILE * f;
  char fname[60];
  strcpy(fname, dir);
  strcat(fname, "/Constants/constants");
  strcat(fname, app);
  strcat(fname, ".dat\0");  f = fopen(fname, "w");

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


void SaveData::saveVar(string variable, int num)
{
  int cpa(0); // cons=1,prims=2,aux=3
  int Nvar(0); // Variable number
  FILE * f;
  char fname[60];

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
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int k(0); k < d->Nz; k++) {
        if (cpa==1) fprintf(f, "%.16f ", d->cons[ID(Nvar, i, j, k)]);
        else if (cpa==2) fprintf(f, "%.16f ", d->prims[ID(Nvar, i, j, k)]);
        else if (cpa==3) fprintf(f, "%.16f ", d->aux[ID(Nvar, i, j, k)]);
        else {
          printf("Error: cpa not set correctly\n");
          exit(1);
        }
      }
      fprintf(f, "\n");
    }
  }

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
