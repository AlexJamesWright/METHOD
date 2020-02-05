#ifndef HYBRID_H
#define HYBRID_H

#include "model.h"
#include "srrmhd.h"
#include "srmhd.h"
#include "REGIME.h"
#include "flux.h"

class Hybrid : public Model
{

  public:

    double
    *icons, *iprims, *iaux,                     // For every cell
    *sicons, *siprims, *siaux,                  // For single cells
    *iflux, *rflux,                             // For every cell
    *isource, *rsource, *regimeSource,          // For every cell
    sigmaCrossOver,
    sigmaSpan;

    bool useREGIME;

    SRRMHD * resistiveModel;

    SRMHD * idealModel;

    REGIME * subgridModel = NULL;


    Hybrid(); //!< Default constructor


    Hybrid(Data * data, double sigmaCrossOver=400, double sigmaSpan=350, bool useREGIME=true);


    ~Hybrid();  //!< Destructor


    void setSubgridModel(FluxMethod * fluxMethod);

    // cons prims aux are all for a single cell!
    double idealWeight(double * cons, double * prims, double * aux);

    // cons prims aux are for all cells!
    double idealWeightID(double * cons, double * prims, double * aux, int i, int j, int k);

    // cons prims aux are all for a single cell!
    bool useResistive(double * cons, double * prims, double * aux);

    // cons prims aux are all for a single cell!
    void setIdealCPAs(double *cons, double * prims, double * aux);

    // cons prims aux are for all cells!
    void setIdealCPAsAll(double *cons, double * prims, double * aux);


    void sourceTermSingleCell(double *cons, double *prims, double *aux, double *source, int i=-1, int j=-1, int k=-1);


    void sourceTerm(double *cons, double *prims, double *aux, double *source);


    void getPrimitiveVarsSingleCell(double *cons, double *prims, double *aux, int i=-1, int j=-1, int k=-1);


    void getPrimitiveVars(double *cons, double *prims, double *aux);


    void primsToAll(double *cons, double *prims, double *aux);


    void fluxVector(double *cons, double *prims, double *aux, double *f, const int dir);


    void finalise(double *cons, double *prims, double *aux);
};

#endif
