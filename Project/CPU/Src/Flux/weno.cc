#include "weno.h"
#include "wenoUpwinds.h"

#define IDD(variable, idx, jdx, kdx) ((variable)*(data->Nx)*(data->Ny)*(data->Nz) + (idx)*(data->Ny)*(data->Nz) + (jdx)*(data->Nz) + (kdx))


WenoBase::WenoBase(Data * data, int order) : data(data), order(order)
{
  shift = (order+1)/2;
  checkSufficientGhostZones();
}

void WenoBase::checkSufficientGhostZones()
{
  if (data->Ng < this->shift+1)
  {
    printf("This order Weno reconstruction requires at least %d ghost zones, you have %d...", shift+1, data->Ng);
    throw std::invalid_argument("You must increase number of boundary cells.");
  }
}


void WenoBase::reconstructUpwind(double * arr, double * recon, int nvars, int dir)
{
  // Syntax
  Data * d(this->data);

  int shft(shift);

  // Reconstruct to determine the flux at the cell face and compute difference
  if (dir == 0) { // x-direction
    # pragma omp parallel for   default(none)   shared(recon, shft, arr, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=shft; i < d->Nx-shft; i++) {
        for (int j=0; j < d->Ny; j++) {
          for (int k=0; k < d->Nz; k++) {
            recon[ID(var, i, j, k)] = upwindX(arr, var, i, j, k);
          }
        }
      }
    }
  }
  else if (dir == 1) { // y-direction
    # pragma omp parallel for   default(none)   shared(recon, arr, shft, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=0; i < d->Nx; i++) {
        for (int j=shft; j < d->Ny-shft; j++) {
          for (int k=0; k < d->Nz; k++) {
            recon[ID(var, i, j, k)] = upwindY(arr, var, i, j, k);
          }
        }
      }
    }
  }
  else { // z-direction
    # pragma omp parallel for   default(none)   shared(recon, arr, shft, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=0; i < d->Nx; i++) {
        for (int j=0; j < d->Ny; j++) {
          for (int k=shft; k < d->Nz-shft; k++) {
            recon[ID(var, i, j, k)] = upwindZ(arr, var, i, j, k);
          }
        }
      }
    }
  }
}

void WenoBase::reconstructDownwind(double * arr, double * recon, int nvars, int dir)
{
  // Syntax
  Data * d(this->data);

  int shft(shift);

  // Reconstruct to determine the flux at the cell face and compute difference
  if (dir == 0) { // x-direction
    # pragma omp parallel for   default (none)   shared (recon, shft, arr, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=shft; i < d->Nx-shft; i++) {
        for (int j=0; j < d->Ny; j++) {
          for (int k=0; k < d->Nz; k++) {
            recon[ID(var, i, j, k)] = downwindX(arr, var, i, j, k);
          }
        }
      }
    }
  }
  else if (dir == 1) { // y-direction
    # pragma omp parallel for   default (none)   shared (recon, arr, shft, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=0; i < d->Nx; i++) {
        for (int j=shft; j < d->Ny-shft; j++) {
          for (int k=0; k < d->Nz; k++) {
            recon[ID(var, i, j, k)] = downwindY(arr, var, i, j, k);
          }
        }
      }
    }
  }
  else { // z-direction
    # pragma omp parallel for   default (none)   shared (recon, arr, shft, d, nvars)
    for (int var=0; var < nvars; var++) {
      for (int i=0; i < d->Nx; i++) {
        for (int j=0; j < d->Ny; j++) {
          for (int k=shft; k < d->Nz-shft; k++) {
            recon[ID(var, i, j, k)] = downwindZ(arr, var, i, j, k);
          }
        }
      }
    }
  }
}




////////////////////////////////////////////////////////////////////////////////
// Weno3
////////////////////////////////////////////////////////////////////////////////


double Weno3::upwindX(double * arr, int var, int i, int j, int k)
{
  return weno3_upwind(arr[IDD(var, i-shift, j, k)],
                      arr[IDD(var, i-shift+1, j, k)],
                      arr[IDD(var, i-shift+2, j, k)]);
}

double Weno3::upwindY(double * arr, int var, int i, int j, int k)
{
  return weno3_upwind(arr[IDD(var, i, j-shift, k)],
                      arr[IDD(var, i, j-shift+1, k)],
                      arr[IDD(var, i, j-shift+2, k)]);
}

double Weno3::upwindZ(double * arr, int var, int i, int j, int k)
{
  return weno3_upwind(arr[IDD(var, i, j, k-shift)],
                      arr[IDD(var, i, j, k-shift+1)],
                      arr[IDD(var, i, j, k-shift+2)]);
}


double Weno3::downwindX(double * arr, int var, int i, int j, int k)
{
  return weno3_upwind(arr[IDD(var, i+shift-1, j, k)],
                      arr[IDD(var, i+shift-2, j, k)],
                      arr[IDD(var, i+shift-3, j, k)]);
}

double Weno3::downwindY(double * arr, int var, int i, int j, int k)
{
  return weno3_upwind(arr[IDD(var, i, j+shift-1, k)],
                      arr[IDD(var, i, j+shift-2, k)],
                      arr[IDD(var, i, j+shift-3, k)]);
}

double Weno3::downwindZ(double * arr, int var, int i, int j, int k)
{
  return weno3_upwind(arr[IDD(var, i, j, k+shift-1)],
                      arr[IDD(var, i, j, k+shift-2)],
                      arr[IDD(var, i, j, k+shift-3)]);
}





////////////////////////////////////////////////////////////////////////////////
// Weno5
////////////////////////////////////////////////////////////////////////////////


double Weno5::upwindX(double * arr, int var, int i, int j, int k)
{
  return weno5_upwind(arr[IDD(var, i-shift  , j, k)],
                      arr[IDD(var, i-shift+1, j, k)],
                      arr[IDD(var, i-shift+2, j, k)],
                      arr[IDD(var, i-shift+3, j, k)],
                      arr[IDD(var, i-shift+4, j, k)]);
}

double Weno5::upwindY(double * arr, int var, int i, int j, int k)
{
  return weno5_upwind(arr[IDD(var, i, j-shift, k)],
                      arr[IDD(var, i, j-shift+1, k)],
                      arr[IDD(var, i, j-shift+2, k)],
                      arr[IDD(var, i, j-shift+3, k)],
                      arr[IDD(var, i, j-shift+4, k)]);
}

double Weno5::upwindZ(double * arr, int var, int i, int j, int k)
{
  return weno5_upwind(arr[IDD(var, i, j, k-shift)],
                      arr[IDD(var, i, j, k-shift+1)],
                      arr[IDD(var, i, j, k-shift+2)],
                      arr[IDD(var, i, j, k-shift+3)],
                      arr[IDD(var, i, j, k-shift+4)]);
}


double Weno5::downwindX(double * arr, int var, int i, int j, int k)
{
  return weno5_upwind(arr[IDD(var, i+shift-1, j, k)],
                      arr[IDD(var, i+shift-2, j, k)],
                      arr[IDD(var, i+shift-3, j, k)],
                      arr[IDD(var, i+shift-4, j, k)],
                      arr[IDD(var, i+shift-5, j, k)]);
}

double Weno5::downwindY(double * arr, int var, int i, int j, int k)
{
  return weno5_upwind(arr[IDD(var, i, j+shift-1, k)],
                      arr[IDD(var, i, j+shift-2, k)],
                      arr[IDD(var, i, j+shift-3, k)],
                      arr[IDD(var, i, j+shift-4, k)],
                      arr[IDD(var, i, j+shift-5, k)]);
}

double Weno5::downwindZ(double * arr, int var, int i, int j, int k)
{
  return weno5_upwind(arr[IDD(var, i, j, k+shift-1)],
                      arr[IDD(var, i, j, k+shift-2)],
                      arr[IDD(var, i, j, k+shift-3)],
                      arr[IDD(var, i, j, k+shift-4)],
                      arr[IDD(var, i, j, k+shift-5)]);
}





////////////////////////////////////////////////////////////////////////////////
// Weno7
////////////////////////////////////////////////////////////////////////////////


double Weno7::upwindX(double * arr, int var, int i, int j, int k)
{
  return weno7_upwind(arr[IDD(var, i-shift  , j, k)],
                      arr[IDD(var, i-shift+1, j, k)],
                      arr[IDD(var, i-shift+2, j, k)],
                      arr[IDD(var, i-shift+3, j, k)],
                      arr[IDD(var, i-shift+4, j, k)],
                      arr[IDD(var, i-shift+5, j, k)],
                      arr[IDD(var, i-shift+6, j, k)]);
}

double Weno7::upwindY(double * arr, int var, int i, int j, int k)
{
  return weno7_upwind(arr[IDD(var, i, j-shift, k)],
                      arr[IDD(var, i, j-shift+1, k)],
                      arr[IDD(var, i, j-shift+2, k)],
                      arr[IDD(var, i, j-shift+3, k)],
                      arr[IDD(var, i, j-shift+4, k)],
                      arr[IDD(var, i, j-shift+5, k)],
                      arr[IDD(var, i, j-shift+6, k)]);
}

double Weno7::upwindZ(double * arr, int var, int i, int j, int k)
{
  return weno7_upwind(arr[IDD(var, i, j, k-shift)],
                      arr[IDD(var, i, j, k-shift+1)],
                      arr[IDD(var, i, j, k-shift+2)],
                      arr[IDD(var, i, j, k-shift+3)],
                      arr[IDD(var, i, j, k-shift+4)],
                      arr[IDD(var, i, j, k-shift+5)],
                      arr[IDD(var, i, j, k-shift+6)]);
}


double Weno7::downwindX(double * arr, int var, int i, int j, int k)
{
  return weno7_upwind(arr[IDD(var, i+shift-1, j, k)],
                      arr[IDD(var, i+shift-2, j, k)],
                      arr[IDD(var, i+shift-3, j, k)],
                      arr[IDD(var, i+shift-4, j, k)],
                      arr[IDD(var, i+shift-5, j, k)],
                      arr[IDD(var, i+shift-6, j, k)],
                      arr[IDD(var, i+shift-7, j, k)]);
}

double Weno7::downwindY(double * arr, int var, int i, int j, int k)
{
  return weno7_upwind(arr[IDD(var, i, j+shift-1, k)],
                      arr[IDD(var, i, j+shift-2, k)],
                      arr[IDD(var, i, j+shift-3, k)],
                      arr[IDD(var, i, j+shift-4, k)],
                      arr[IDD(var, i, j+shift-5, k)],
                      arr[IDD(var, i, j+shift-6, k)],
                      arr[IDD(var, i, j+shift-7, k)]);
}

double Weno7::downwindZ(double * arr, int var, int i, int j, int k)
{
  return weno7_upwind(arr[IDD(var, i, j, k+shift-1)],
                      arr[IDD(var, i, j, k+shift-2)],
                      arr[IDD(var, i, j, k+shift-3)],
                      arr[IDD(var, i, j, k+shift-4)],
                      arr[IDD(var, i, j, k+shift-5)],
                      arr[IDD(var, i, j, k+shift-6)],
                      arr[IDD(var, i, j, k+shift-7)]);
}




////////////////////////////////////////////////////////////////////////////////
// Weno9
////////////////////////////////////////////////////////////////////////////////


double Weno9::upwindX(double * arr, int var, int i, int j, int k)
{
  return weno9_upwind(arr[IDD(var, i-shift  , j, k)],
                      arr[IDD(var, i-shift+1, j, k)],
                      arr[IDD(var, i-shift+2, j, k)],
                      arr[IDD(var, i-shift+3, j, k)],
                      arr[IDD(var, i-shift+4, j, k)],
                      arr[IDD(var, i-shift+5, j, k)],
                      arr[IDD(var, i-shift+6, j, k)],
                      arr[IDD(var, i-shift+7, j, k)],
                      arr[IDD(var, i-shift+8, j, k)]);
}

double Weno9::upwindY(double * arr, int var, int i, int j, int k)
{
  return weno9_upwind(arr[IDD(var, i, j-shift, k)],
                      arr[IDD(var, i, j-shift+1, k)],
                      arr[IDD(var, i, j-shift+2, k)],
                      arr[IDD(var, i, j-shift+3, k)],
                      arr[IDD(var, i, j-shift+4, k)],
                      arr[IDD(var, i, j-shift+5, k)],
                      arr[IDD(var, i, j-shift+6, k)],
                      arr[IDD(var, i, j-shift+7, k)],
                      arr[IDD(var, i, j-shift+8, k)]);
}

double Weno9::upwindZ(double * arr, int var, int i, int j, int k)
{
  return weno9_upwind(arr[IDD(var, i, j, k-shift)],
                      arr[IDD(var, i, j, k-shift+1)],
                      arr[IDD(var, i, j, k-shift+2)],
                      arr[IDD(var, i, j, k-shift+3)],
                      arr[IDD(var, i, j, k-shift+4)],
                      arr[IDD(var, i, j, k-shift+5)],
                      arr[IDD(var, i, j, k-shift+6)],
                      arr[IDD(var, i, j, k-shift+7)],
                      arr[IDD(var, i, j, k-shift+8)]);
}


double Weno9::downwindX(double * arr, int var, int i, int j, int k)
{
  return weno9_upwind(arr[IDD(var, i+shift-1, j, k)],
                      arr[IDD(var, i+shift-2, j, k)],
                      arr[IDD(var, i+shift-3, j, k)],
                      arr[IDD(var, i+shift-4, j, k)],
                      arr[IDD(var, i+shift-5, j, k)],
                      arr[IDD(var, i+shift-6, j, k)],
                      arr[IDD(var, i+shift-7, j, k)],
                      arr[IDD(var, i+shift-8, j, k)],
                      arr[IDD(var, i+shift-9, j, k)]);
}

double Weno9::downwindY(double * arr, int var, int i, int j, int k)
{
  return weno9_upwind(arr[IDD(var, i, j+shift-1, k)],
                      arr[IDD(var, i, j+shift-2, k)],
                      arr[IDD(var, i, j+shift-3, k)],
                      arr[IDD(var, i, j+shift-4, k)],
                      arr[IDD(var, i, j+shift-5, k)],
                      arr[IDD(var, i, j+shift-6, k)],
                      arr[IDD(var, i, j+shift-7, k)],
                      arr[IDD(var, i, j+shift-8, k)],
                      arr[IDD(var, i, j+shift-9, k)]);
}

double Weno9::downwindZ(double * arr, int var, int i, int j, int k)
{
  return weno9_upwind(arr[IDD(var, i, j, k+shift-1)],
                      arr[IDD(var, i, j, k+shift-2)],
                      arr[IDD(var, i, j, k+shift-3)],
                      arr[IDD(var, i, j, k+shift-4)],
                      arr[IDD(var, i, j, k+shift-5)],
                      arr[IDD(var, i, j, k+shift-6)],
                      arr[IDD(var, i, j, k+shift-7)],
                      arr[IDD(var, i, j, k+shift-8)],
                      arr[IDD(var, i, j, k+shift-9)]);
}






////////////////////////////////////////////////////////////////////////////////
// Weno11
////////////////////////////////////////////////////////////////////////////////


double Weno11::upwindX(double * arr, int var, int i, int j, int k)
{
  return weno11_upwind(arr[IDD(var, i-shift  , j, k)],
                      arr[IDD(var, i-shift+1, j, k)],
                      arr[IDD(var, i-shift+2, j, k)],
                      arr[IDD(var, i-shift+3, j, k)],
                      arr[IDD(var, i-shift+4, j, k)],
                      arr[IDD(var, i-shift+5, j, k)],
                      arr[IDD(var, i-shift+6, j, k)],
                      arr[IDD(var, i-shift+7, j, k)],
                      arr[IDD(var, i-shift+8, j, k)],
                      arr[IDD(var, i-shift+9, j, k)],
                      arr[IDD(var, i-shift+10, j, k)]);
}

double Weno11::upwindY(double * arr, int var, int i, int j, int k)
{
  return weno11_upwind(arr[IDD(var, i, j-shift, k)],
                      arr[IDD(var, i, j-shift+1, k)],
                      arr[IDD(var, i, j-shift+2, k)],
                      arr[IDD(var, i, j-shift+3, k)],
                      arr[IDD(var, i, j-shift+4, k)],
                      arr[IDD(var, i, j-shift+5, k)],
                      arr[IDD(var, i, j-shift+6, k)],
                      arr[IDD(var, i, j-shift+7, k)],
                      arr[IDD(var, i, j-shift+8, k)],
                      arr[IDD(var, i, j-shift+9, k)],
                      arr[IDD(var, i, j-shift+10, k)]);
}

double Weno11::upwindZ(double * arr, int var, int i, int j, int k)
{
  return weno11_upwind(arr[IDD(var, i, j, k-shift)],
                      arr[IDD(var, i, j, k-shift+1)],
                      arr[IDD(var, i, j, k-shift+2)],
                      arr[IDD(var, i, j, k-shift+3)],
                      arr[IDD(var, i, j, k-shift+4)],
                      arr[IDD(var, i, j, k-shift+5)],
                      arr[IDD(var, i, j, k-shift+6)],
                      arr[IDD(var, i, j, k-shift+7)],
                      arr[IDD(var, i, j, k-shift+8)],
                      arr[IDD(var, i, j, k-shift+9)],
                      arr[IDD(var, i, j, k-shift+10)]);
}


double Weno11::downwindX(double * arr, int var, int i, int j, int k)
{
  return weno11_upwind(arr[IDD(var, i+shift-1, j, k)],
                      arr[IDD(var, i+shift-2, j, k)],
                      arr[IDD(var, i+shift-3, j, k)],
                      arr[IDD(var, i+shift-4, j, k)],
                      arr[IDD(var, i+shift-5, j, k)],
                      arr[IDD(var, i+shift-6, j, k)],
                      arr[IDD(var, i+shift-7, j, k)],
                      arr[IDD(var, i+shift-8, j, k)],
                      arr[IDD(var, i+shift-9, j, k)],
                      arr[IDD(var, i+shift-10, j, k)],
                      arr[IDD(var, i+shift-11, j, k)]);
}

double Weno11::downwindY(double * arr, int var, int i, int j, int k)
{
  return weno11_upwind(arr[IDD(var, i, j+shift-1, k)],
                      arr[IDD(var, i, j+shift-2, k)],
                      arr[IDD(var, i, j+shift-3, k)],
                      arr[IDD(var, i, j+shift-4, k)],
                      arr[IDD(var, i, j+shift-5, k)],
                      arr[IDD(var, i, j+shift-6, k)],
                      arr[IDD(var, i, j+shift-7, k)],
                      arr[IDD(var, i, j+shift-8, k)],
                      arr[IDD(var, i, j+shift-9, k)],
                      arr[IDD(var, i, j+shift-10, k)],
                      arr[IDD(var, i, j+shift-11, k)]);
}

double Weno11::downwindZ(double * arr, int var, int i, int j, int k)
{
  return weno11_upwind(arr[IDD(var, i, j, k+shift-1)],
                      arr[IDD(var, i, j, k+shift-2)],
                      arr[IDD(var, i, j, k+shift-3)],
                      arr[IDD(var, i, j, k+shift-4)],
                      arr[IDD(var, i, j, k+shift-5)],
                      arr[IDD(var, i, j, k+shift-6)],
                      arr[IDD(var, i, j, k+shift-7)],
                      arr[IDD(var, i, j, k+shift-8)],
                      arr[IDD(var, i, j, k+shift-9)],
                      arr[IDD(var, i, j, k+shift-10)],
                      arr[IDD(var, i, j, k+shift-11)]);
}
