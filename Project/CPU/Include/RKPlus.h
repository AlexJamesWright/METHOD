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


//! Second order RK
/*!
    SSPRK(2,2) from Gottlieb 2009
    This can run the AdvectionSingleFluid test with cfl=0.6
*/
class RK2B : public RKPlus
{
  public:

    // Need some work arrays
    double *u1cons, *u1prims, *u1aux, *rhs1, *rhs2;

    RK2B(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    ~RK2B();

    virtual void stage1(double * cons, double * prims, double * aux, double dt);

    void stage2(double * cons, double * prims, double * aux, double dt);

    void step(double * cons, double * prims, double * aux, double dt=0);

};

//! Third order RK
/*!
    SSPRK(3,3) from Gottlieb 2009
    This can run the AdvectionSingleFluid test with cfl=1.7
*/
class RK3 : public RKPlus
{
  public:

    // Need some work arrays
    double *u1cons, *u1prims, *u1aux, *u2cons, *u2prims, *u2aux, *rhs1, *rhs2, *rhs3;

    RK3(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    ~RK3();

    void stage1(double * cons, double * prims, double * aux, double dt);

    void stage2(double * cons, double * prims, double * aux, double dt);

    void stage3(double * cons, double * prims, double * aux, double dt);

    void step(double * cons, double * prims, double * aux, double dt=0);

};

//! Fourth order RK
/*!
    SSPRK(5,4) from Gottlieb 2009
    This can run the AdvectionSingleFluid test with cfl=2.3
*/
class RK4 : public RKPlus
{
  public:

    // Need some work arrays
    double *u1cons, *u1prims, *u1aux, *u2cons, *u2prims, *u2aux, *u3cons, *u3prims, *u3aux, *u4cons, *u4prims, *u4aux, *rhs1, *rhs2, *rhs3, *rhs4, *rhs5;

    RK4(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    ~RK4();

    void stage1(double * cons, double * prims, double * aux, double dt);

    void stage2(double * cons, double * prims, double * aux, double dt);

    void stage3(double * cons, double * prims, double * aux, double dt);

    void stage4(double * cons, double * prims, double * aux, double dt);

    void stage5(double * cons, double * prims, double * aux, double dt);

    void step(double * cons, double * prims, double * aux, double dt=0);

};


//! Fourth order RK
/*!
    SSPRK(10,4) from Gottlieb 2009
    This can run the AdvectionSingleFluid test with cfl=3.6
*/
class RK4_10 : public RKPlus
{
  public:

    // Need some work arrays
    double *u1cons, *u1prims, *u1aux, *u2cons, *u2prims, *u2aux, *rhs1;

    RK4_10(Data * data, Model * model, Bcs * bcs, FluxMethod * fluxMethod, ModelExtension * modelExtension = NULL);

    ~RK4_10();

    void prepare1(double * cons, double * prims, double * aux);

    void prepare2(double * cons, double * prims, double * aux);

    void stageRepeat(double * cons, double * prims, double * aux, double dt);

    void stageFinal(double * cons, double * prims, double * aux, double dt);

    void step(double * cons, double * prims, double * aux, double dt=0);

};


#endif
