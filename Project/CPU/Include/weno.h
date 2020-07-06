#ifndef WENO_H
#define WENO_H

#include "simData.h"
#include <stdexcept>

class WenoBase
{
  public:

    Data * data;

    int order;    //!< Order of reconstruction and number of buffers for this scheme (should be 1 more than data->Ng)
    int shift;    //!< Shift = (order+1)/2

    WenoBase(Data * data) : data(data) { }

    virtual void reconstructUpwind(double * arr, double * recon, int nvars, int dir) = 0;

    virtual void reconstructDownwind(double * arr, double * recon, int nvars, int dir) = 0;

    virtual void checkSufficientGhostZones();

};


class Weno3 : public WenoBase
{
  public:
    Weno3(Data * data) : WenoBase(data)
    {
      this->order = 3;
      this->shift = (order+1)/2;
      checkSufficientGhostZones();
    }

    virtual void reconstructUpwind(double * arr, double * recon, int nvars, int dir);

    virtual void reconstructDownwind(double * arr, double * recon, int nvars, int dir);

};


#endif
