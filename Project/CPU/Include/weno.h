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

    WenoBase(Data * data, int order);

    virtual double upwindX(double * arr, int var, int i, int j, int k) = 0;
    virtual double upwindY(double * arr, int var, int i, int j, int k) = 0;
    virtual double upwindZ(double * arr, int var, int i, int j, int k) = 0;

    virtual double downwindX(double * arr, int var, int i, int j, int k) = 0;
    virtual double downwindY(double * arr, int var, int i, int j, int k) = 0;
    virtual double downwindZ(double * arr, int var, int i, int j, int k) = 0;

    virtual void reconstructUpwind(double * arr, double * recon, int nvars, int dir);

    virtual void reconstructDownwind(double * arr, double * recon, int nvars, int dir);

    virtual void checkSufficientGhostZones();

};


class Weno3 : public WenoBase
{
  public:
    Weno3(Data * data) : WenoBase(data, 3) { }

    virtual double upwindX(double * arr, int var, int i, int j, int k);
    virtual double upwindY(double * arr, int var, int i, int j, int k);
    virtual double upwindZ(double * arr, int var, int i, int j, int k);

    virtual double downwindX(double * arr, int var, int i, int j, int k);
    virtual double downwindY(double * arr, int var, int i, int j, int k);
    virtual double downwindZ(double * arr, int var, int i, int j, int k);

};


class Weno5 : public WenoBase
{
  public:
    Weno5(Data * data) : WenoBase(data, 5) { }

    virtual double upwindX(double * arr, int var, int i, int j, int k);
    virtual double upwindY(double * arr, int var, int i, int j, int k);
    virtual double upwindZ(double * arr, int var, int i, int j, int k);

    virtual double downwindX(double * arr, int var, int i, int j, int k);
    virtual double downwindY(double * arr, int var, int i, int j, int k);
    virtual double downwindZ(double * arr, int var, int i, int j, int k);

};

class Weno7 : public WenoBase
{
  public:
    Weno7(Data * data) : WenoBase(data, 7) { }

    virtual double upwindX(double * arr, int var, int i, int j, int k);
    virtual double upwindY(double * arr, int var, int i, int j, int k);
    virtual double upwindZ(double * arr, int var, int i, int j, int k);

    virtual double downwindX(double * arr, int var, int i, int j, int k);
    virtual double downwindY(double * arr, int var, int i, int j, int k);
    virtual double downwindZ(double * arr, int var, int i, int j, int k);

};

class Weno9 : public WenoBase
{
  public:
    Weno9(Data * data) : WenoBase(data, 9) { }

    virtual double upwindX(double * arr, int var, int i, int j, int k);
    virtual double upwindY(double * arr, int var, int i, int j, int k);
    virtual double upwindZ(double * arr, int var, int i, int j, int k);

    virtual double downwindX(double * arr, int var, int i, int j, int k);
    virtual double downwindY(double * arr, int var, int i, int j, int k);
    virtual double downwindZ(double * arr, int var, int i, int j, int k);

};

class Weno11 : public WenoBase
{
  public:
    Weno11(Data * data) : WenoBase(data, 11) { }

    virtual double upwindX(double * arr, int var, int i, int j, int k);
    virtual double upwindY(double * arr, int var, int i, int j, int k);
    virtual double upwindZ(double * arr, int var, int i, int j, int k);

    virtual double downwindX(double * arr, int var, int i, int j, int k);
    virtual double downwindY(double * arr, int var, int i, int j, int k);
    virtual double downwindZ(double * arr, int var, int i, int j, int k);

};

#endif
