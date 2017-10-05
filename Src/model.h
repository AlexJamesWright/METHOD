#ifndef MODEL_H
#define MODEL_H

#include "simulation.h"

/*
    We're using an abstract base class to form the foundations of the models preparing
  for additional set ups. Systems will derive from Model and supply their own
  spectral analyses for cons2prims, prims2cons functions and the flux vectors.
  Constructors must require access to public Simulation data.
*/

class Model
{
  protected:
    Simulation *sim;  // Pointer to overarching Simulation object
    int Ncons, Nprims, Naux;  // Size of conserved, primitive and aux state vectors

  public:
    Model();
    Model(Simulation *sim) : sim(sim) {}
    ~Model() {}
    virtual void fluxFunc(double *cons, double *prims, double *aux, int dir) = 0;
    virtual void getPrimitiveVars(double *cons, double *prims, double *aux) = 0;
    virtual void consToPrims(double *cons, double *prims, double *aux) = 0;

};

class SRMHD : public Model
{
  public:
    // Constructor
    SRMHD();
    SRMHD(Simulation *sim);
    ~SRMHD();

    // Member functions
    void fluxFunc(double *cons, double *prims, double *aux, int dir);
    void getPrimitiveVars(double *cons, double *prims, double *aux);
    void consToPrims(double *cons, double *prims, double *aux);
};

#endif
