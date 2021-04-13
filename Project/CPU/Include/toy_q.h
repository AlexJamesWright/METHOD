#ifndef TOY_Q_H
#define TOY_Q_H

#include "model.h"
#include "simData.h"

//! Toy model for (relativistic) heat equation
/*!
  @par
    Toy model for heat equation: essentially, this is Cattaneo.
*/
class ToyQ : public Model
{
  public:
    ToyQ();
    ToyQ(Data * data);
    virtual ~ToyQ();

    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);
    void sourceTerm(double *cons, double *prims, double *aux, double *source);
    void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1);
    void getPrimitiveVars(double *cons, double *prims, double *aux);
    void primsToAll(double *cons, double *prims, double *aux);
    void fluxVector(double *cons, double *prims, double *aux, double *f, const int dir);
    void finalise(double *cons, double *prims, double *aux) { };



};

//! Toy model for (relativistic) heat equation
/*!
  @par
    Toy model for heat equation: essentially, this is Cattaneo.

    Includes (hard coded) functional dependence of kappa, tau on T.
*/
class ToyQFunctional : public ToyQ
{
  public:
    ToyQFunctional();
    ToyQFunctional(Data * data);
    virtual ~ToyQFunctional();

    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);
    void sourceTerm(double *cons, double *prims, double *aux, double *source);

};

#endif
