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

void ResistiveSGM::reset(void)
{
  // TODO
}

void ResistiveSGM::set_vars(double * cons, double * prims, double * aux)
{
  // TODO
}

void ResistiveSGM::set_dwdsb(double * cons, double * prims, double * aux)
{
  // TODO
}

void ResistiveSGM::set_Dx(double * cons, double * prims, double * aux)
{
  // TODO
}

void ResistiveSGM::set_Dy(double * cons, double * prims, double * aux)
{
  // TODO
}


void ResistiveSGM::set_Dz(double * cons, double * prims, double * aux)
{
  // TODO
}


void ResistiveSGM::set_K(double * cons, double * prims, double * aux)
{
  // TODO
}


void ResistiveSGM::set_dfxdw(double * cons, double * prims, double * aux)
{
  // TODO
}


void ResistiveSGM::set_dfydw(double * cons, double * prims, double * aux)
{
  // TODO
}


void ResistiveSGM::set_dfzdw(double * cons, double * prims, double * aux)
{
  // TODO
}
