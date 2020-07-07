#ifndef RKPLUS_H
#define RKPLUS_H

#include "timeInt.h"
/*!
    For the higher order RK intergrators, we are going to bundle the source
  and the flux into a RHS function
*/
class RKPlus : public TimeIntegrator
{
  public:

    double * fluxCont;

    RKPlus(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    ~RKPlus();

    virtual void rhs(double * cons, double * prims, double * aux, double * rhsVec);

};


//! RK2B integrator, but based off RKPlus API
class RK2B : public RKPlus
{
  public:
    // Need some work arrays
    double *p1cons, *p1prims, *p1aux, *args1, *args2;


    RK2B(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    ~RK2B();

    void predictorStep(double * cons, double * prims, double * aux, double dt);

    void correctorStep(double * cons, double * prims, double * aux, double dt);

    void step(double * cons, double * prims, double * aux, double dt=0);


};

#endif
