#include "model.h"

SRMHD::SRMHD() : Model()
{
  this->Ncons = 9;
  this->Nprims = 8;
  this->Naux = 9;
}

SRMHD::SRMHD(Simulation * sim) : Model(sim)
{
  this->Ncons = (this->sim)->Ncons = 9;
  this->Nprims = (this->sim)->Nprims = 8;
  this->Naux = (this->sim)->Naux = 9;
}


void SRMHD::fluxFunc(double *cons, double *prims, double *aux, int dir)
{

}

void SRMHD::getPrimitiveVars(double *cons, double *prims, double *aux)
{

}

void SRMHD::consToPrims(double *cons, double *prims, double *aux)
{

}
