#ifndef TOY_Q_H
#define TOY_Q_H

#include "model.h"
#include "simData.h"

//! Non relativistic Euler equations
/*!
  @par
    Standard, idea fluid, non relativistic euler equations, taken from Leveque 14.3
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

#endif
