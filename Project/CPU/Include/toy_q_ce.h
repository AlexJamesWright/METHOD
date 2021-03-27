#ifndef TOY_Q_CE_H
#define TOY_Q_CE_H

#include "model.h"
#include "simData.h"

//! Toy model for (relativistic) heat equation in Chapman Enskog limit
/*!
  @par
    Toy model for heat equation: essentially, this is Cattaneo.
*/
class ToyQ_CE : public Model
{
  public:
    ToyQ_CE();
    ToyQ_CE(Data * data);
    virtual ~ToyQ_CE();

    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);
    void sourceTerm(double *cons, double *prims, double *aux, double *source);
    void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1);
    void getPrimitiveVars(double *cons, double *prims, double *aux);
    void primsToAll(double *cons, double *prims, double *aux);
    void fluxVector(double *cons, double *prims, double *aux, double *f, const int dir);
    void finalise(double *cons, double *prims, double *aux) { };



};

#endif
