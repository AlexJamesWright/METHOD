#include "initFunc.h"
#include <iostream>

InitialFunc::InitialFunc(Data * data) : data(data)
{
  Data * d;
  d = this->data;

  std::cout << "Element (4, 6, 6) = " << d->cons[4*d->Ny*d->Nx + 6*d->Ny + 6] << std::endl;

  std::cout << "Setting cons data..." << std::endl;

  // Set all state vectors to zero
  for (int i(0); i < d->Nx; i++) {
    for (int j(0); j < d->Ny; j++) {
      for (int var(0); var < d->Ncons; var++) {
        d->cons[(var, i, j)] = j;
      }
    }
  }

  std::cout << "Element (4, 6, 6) = " << d->cons[4*d->Ny*d->Nx + 6*d->Ny + 6] << std::endl;

}

OTVortex::OTVortex(Data * data) : InitialFunc(data)
{

}
