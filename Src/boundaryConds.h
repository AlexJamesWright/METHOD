#ifndef BOUNDARYCONDS_H
#define BOUNDARYCONDS_H

#include "simulation.h"

  class Bcs
  {

    public:
      Bcs(*Simulation sim);
      void outflow(**double);
      void periodic(**double);

  };

#endif
