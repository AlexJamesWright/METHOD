#include "simData.h"
#include "platformEnv.h"
#include <stdexcept>
#include <cmath>
#include <string>

#define IDn(variable, idx, jdx, kdx) ((variable)*(this->Nx)*(this->Ny)*(this->Nz) + (idx)*(this->Ny)*(this->Nz) + (jdx)*(this->Nz) + (kdx))

using namespace std;

Data::Data(int nx, int ny, int nz,
           double xmin, double xmax,
           double ymin, double ymax,
           double zmin, double zmax,
           double endTime, PlatformEnv *env,
	         double cfl, int Ng,
           double gamma, double sigma,
           double cp,
           double mu1, double mu2, int frameSkip,
           int reportItersPeriod,
           bool functionalSigma, double gam)
           :
           nx(nx), ny(ny), nz(nz),
           xmin(xmin), xmax(xmax),
           ymin(ymin), ymax(ymax),
           zmin(zmin), zmax(zmax),
           endTime(endTime), cfl(cfl), Ng(Ng),
           gamma(gamma), sigma(sigma),
           memSet(0), bcsSet(0),
           Ncons(0), Nprims(0), Naux(0),
           cp(cp),
           mu1(mu1), mu2(mu2),
           frameSkip(frameSkip),
           reportItersPeriod(reportItersPeriod),
           functionalSigma(functionalSigma), gam(gam), 
	   t(0)
{
	initData(env);
}

Data::Data(DataArgsBase args, PlatformEnv *env)
           :
           nx(args.nx), ny(args.ny), nz(args.nz),
           xmin(args.xmin), xmax(args.xmax),
           ymin(args.ymin), ymax(args.ymax),
           zmin(args.zmin), zmax(args.zmax),
           endTime(args.endTime), cfl(args.cfl), Ng(args.Ng),
           gamma(args.gamma), sigma(args.sigma),
           memSet(0), bcsSet(0),
           Ncons(0), Nprims(0), Naux(0),
           cp(args.cp),
           mu1(args.mu1), mu2(args.mu2),
           frameSkip(args.frameSkip),
           reportItersPeriod(args.reportItersPeriod),
           functionalSigma(args.functionalSigma), gam(args.gam),
           t(args.t) 
{
	initData(env, args.nOptionalSimArgs, args.optionalSimArgs, args.optionalSimArgNames);
}

void Data::initData(PlatformEnv *env, int nOptionalSimArgs, std::vector<double> optionalSimArgs, std::vector<std::string> optionalSimArgNames){
  // TODO -- handle nx not dividing perfectly into nxRanks

  // Set Nx to be nx per MPI process + ghost cells
  this->Nx = nx/env->nxRanks + 2 * Ng;
  this->Ny = ny/env->nyRanks + 2 * Ng;
  this->Nz = nz/env->nzRanks + 2 * Ng;
  this->Ntot = this->Nx * this->Ny * this->Nz;

  //printf("proc %d (%d) initialized with %d nx\n", env->rank, env->xRankId, this->Nx);
  dims = 3;

  // Catch 2D case
  if (nz == 0) {
    this->Nz = 1;
    zmin = -1e20;
    zmax = 1e20;
    dims = 2;
  }
  // Catch 1D case
  if (ny == 0) {
    this->Nz = this->Ny = 1;
    zmin = ymin = -1e20;
    zmax = ymax = 1e20;
    dims = 1;
  }

  // Set some variables that define the interior cells
  is = Ng; ie = Nx-Ng;  // i-start, i-end
  js = Ng; je = Ny-Ng;  // j-start, j-end
  ks = Ng; ke = Nz-Ng;  // k-start, k-end
  if (dims<3) {
    ks = 0; ke = 1;
  }
  if (dims<2) {
    js = 0; je = 1;
  }
  for (int l(1); l < Ng; l++) {
    this->is_minus.push_back(is-l);
    this->ie_plus.push_back(ie+l);
    if (dims > 1) {
      this->js_minus.push_back(js-l);
      this->je_plus.push_back(je+l);
      if (dims > 2) {
        this->ks_minus.push_back(ks-l);
        this->ke_plus.push_back(ke+l);
      }
      else {
        this->ks_minus.push_back(0);
        this->ke_plus.push_back(1);
      }
    }
    else {
      this->js_minus.push_back(0);
      this->je_plus.push_back(1);
      this->ks_minus.push_back(0);
      this->ke_plus.push_back(1);
    }
  }

  // Ensure there is some Resistivity
  if (this->sigma < 0.0) {
    throw std::invalid_argument("Conductivity must be non-negative, sigma >= 0.\n");
  }
  // Ensure charges are correct way round
  if (this->mu1 > 0.0 or this->mu2 < 0.0) {
    throw std::invalid_argument("Species 1 must have negative charge, mu1 < 0, and species 2 must have positive charge, mu2 > 0.\n");
  }

  // Allocate and initialise optional simulation parameters if we have been passed any
  this->nOptionalSimArgs = nOptionalSimArgs;
  this->optionalSimArgs = optionalSimArgs;
  this->optionalSimArgNames = optionalSimArgNames;
}

double Data::sigmaFunc(double * cons, double * prims, double * aux, int i, int j, int k)
{
  if (!functionalSigma)
  {
    return sigma;
  }
  else
  {
    if (i < 0 || j < 0 || k < 0)
      return sigma * pow(cons[0], gam);
    else
      return sigma * pow(cons[IDn(0, i, j, k)], gam);
  }
}

