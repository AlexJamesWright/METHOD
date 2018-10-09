#include "resistiveSGM.h"
#include <cstdio>

ResistiveSGM::ResistiveSGM(Data * data, FluxMethod * fluxMethod) : SubGridModel(data), fluxMethod(fluxMethod)
{
  // Allocate arrays
}

void ResistiveSGM::subgridSource(double * cons, double * prims, double * aux, double * source)
{
  // TODO
}
