#include "resistiveSGM.h"
#include <cstdio>

#define ID(variable, idx, jdx, kdx)  ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

// dwdsb
#define IDWS(ldx, mdx, idx, jdx, kdx)  ((ldx)*(3)*(d->Nx)*(d->Ny)*(d->Nz) + (mdx)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))


ResistiveSGM::ResistiveSGM(Data * data, FluxMethod * fluxMethod) : SubGridModel(data), fluxMethod(fluxMethod)
{
  //  Syntax
  Data * d(this->data);

  // Allocate arrays
  dfxdw = new double[d->Nx*d->Ny*d->Nz*9*12] ();
  dfydw = new double[d->Nx*d->Ny*d->Nz*9*12] ();
  dfzdw = new double[d->Nx*d->Ny*d->Nz*9*12] ();
  dwdsb = new double[d->Nx*d->Ny*d->Nz*12*3] ();
  E = new double[d->Nx*d->Ny*d->Nz*3] ();
  q = new double[d->Nx*d->Ny*d->Nz] ();
  K = new double[d->Nx*d->Ny*d->Nz*3] ();
  Mx = new double[d->Nx*d->Ny*d->Nz*9*3] ();
  My = new double[d->Nx*d->Ny*d->Nz*9*3] ();
  Mz = new double[d->Nx*d->Ny*d->Nz*9*3] ();
  fx = new double[d->Nx*d->Ny*d->Nz*3] ();
  fy = new double[d->Nx*d->Ny*d->Nz*3] ();
  fz = new double[d->Nx*d->Ny*d->Nz*3] ();
  diffuX = new double[d->Nx*d->Ny*d->Nz*d->Ncons] ();
  diffuY = new double[d->Nx*d->Ny*d->Nz*d->Ncons] ();
  diffuZ = new double[d->Nx*d->Ny*d->Nz*d->Ncons] ();
  alpha = new double[d->Nx*d->Ny*d->Nz] ();
}

ResistiveSGM::~ResistiveSGM()
{
  // Clean up your mess
  free(dfxdw);
  free(dfydw);
  free(dfzdw);
  free(dwdsb);
  free(E);
  free(q);
  free(K);
  free(Mx);
  free(My);
  free(Mz);
  free(fx);
  free(fy);
  free(fz);
  free(diffuX);
  free(diffuY);
  free(diffuZ);
  free(alpha);
}

void ResistiveSGM::subgridSource(double * cons, double * prims, double * aux, double * source)
{
  // Zero arrays and set vars
  reset(source);
  set_vars(cons, prims, aux);

  // Ensure K and dwdsb are set
  set_K(cons, prims, aux);
  set_dwdsb(cons, prims, aux);

  // TODO
  // Determine Dx (and Dy and Dz if multi-dimensional domain)
  // Apply MidMod to gradient and add to source vector.
}

void ResistiveSGM::reset(double * source)
{
  // Syntax
  Data * d(this->data);
  // Reset the arrays in which we use the += operator
  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        for (int var(0); var<d->Ncons; var++) {
          source[ID(var, i, j, k)] = 0.0;
        }
        K[ID(0, i, j, k)] = 0.0;
      }
    }
  }
}

void ResistiveSGM::set_vars(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  // First set E = - v cross B
  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        E[ID(0, i, j, k)] = - (prims[ID(2, i, j, k)] * prims[ID(7, i, j, k)] - prims[ID(3, i, j, k)] * prims[ID(6, i, j, k)]);
        E[ID(1, i, j, k)] = - (prims[ID(3, i, j, k)] * prims[ID(5, i, j, k)] - prims[ID(1, i, j, k)] * prims[ID(7, i, j, k)]);
        E[ID(2, i, j, k)] = - (prims[ID(1, i, j, k)] * prims[ID(6, i, j, k)] - prims[ID(2, i, j, k)] * prims[ID(5, i, j, k)]);
      }
    }
  }

  // Charge density, q, is the gradient of the electric field. Central differencing
  if (d->dims == 3) {
    for (int i(2); i<d->Nx-2; i++) {
      for (int j(2); j<d->Ny-2; j++) {
        for (int k(2); k<d->Nz-2; k++) {
          q[ID(0, i, j, k)] = (-E[ID(0, i+2, j, k)] + 8*E[ID(0, i+1, j, k)] - 8*E[ID(0, i-1, j, k)] + E[ID(0, i-2, j, k)]) / (12*d->dx) +
                              (-E[ID(1, i, j+2, k)] + 8*E[ID(1, i, j+1, k)] - 8*E[ID(1, i, j-1, k)] + E[ID(1, i, j-2, k)]) / (12*d->dy) +
                              (-E[ID(2, i, j, k+2)] + 8*E[ID(2, i, j, k+1)] - 8*E[ID(2, i, j, k-1)] + E[ID(2, i, j, k-2)]) / (12*d->dz);
        }
      }
    }
  }
  else if (d->dims == 2) {
    for (int i(2); i<d->Nx-2; i++) {
      for (int j(2); j<d->Ny-2; j++) {
        q[ID(0, i, j, 0)] = (-E[ID(0, i+2, j, 0)] + 8*E[ID(0, i+1, j, 0)] - 8*E[ID(0, i-1, j, 0)] + E[ID(0, i-2, j, 0)]) / (12*d->dx) +
                            (-E[ID(1, i, j+2, 0)] + 8*E[ID(1, i, j+1, 0)] - 8*E[ID(1, i, j-1, 0)] + E[ID(1, i, j-2, 0)]) / (12*d->dy);
      }
    }
  }
  else {
    for (int i(2); i<d->Nx-2; i++) {
      q[ID(0, i, 0, 0)] = (-E[ID(0, i+2, 0, 0)] + 8*E[ID(0, i+1, 0, 0)] - 8*E[ID(0, i-1, 0, 0)] + E[ID(0, i-2, 0, 0)]) / (12*d->dx);
    }
  }


  // Set prefactor for dwdsb
  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        double qsigsq(q[ID(0, i, j, k)]*q[ID(0, i, j, k)] + d->sigma*d->sigma);
        double Bsq(prims[ID(5, i, j, k)]*prims[ID(5, i, j, k)] +
                   prims[ID(6, i, j, k)]*prims[ID(6, i, j, k)] +
                   prims[ID(7, i, j, k)]*prims[ID(7, i, j, k)]);
        alpha[ID(0, i, j, k)] = 1 / (qsigsq*(qsigsq + d->sigma*d->sigma*Bsq));
      }
    }
  }

}


void ResistiveSGM::set_dwdsb(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        // Save some typing
        double qsigsq(q[ID(0, i, j, k)]*q[ID(0, i, j, k)] + d->sigma*d->sigma);
        double sigcu(d->sigma*d->sigma*d->sigma);
        double sig(d->sigma);
        double qch(q[ID(0, i, j, k)]);
        double vdotB(prims[ID(1, i, j, k)]*prims[ID(5, i, j, k)] +
                     prims[ID(2, i, j, k)]*prims[ID(6, i, j, k)] +
                     prims[ID(3, i, j, k)]*prims[ID(7, i, j, k)]);


        // First, do A
        {
          dwdsb[IDWS(1, 0, i, j, k)] = -qch * (qch*qch + sig*sig * (1 + prims[ID(5, i, j, k)]*prims[ID(5, i, j, k)]));
          dwdsb[IDWS(1, 1, i, j, k)] = -sig * (prims[ID(5, i, j, k)] * prims[ID(6, i, j, k)]*qch*sig - prims[ID(7, i, j, k)]*qsigsq);
          dwdsb[IDWS(1, 2, i, j, k)] = -sig * (prims[ID(5, i, j, k)] * prims[ID(7, i, j, k)]*qch*sig + prims[ID(6, i, j, k)]*qsigsq);
          dwdsb[IDWS(2, 0, i, j, k)] = -sig * (prims[ID(5, i, j, k)] * prims[ID(6, i, j, k)]*qch*sig + prims[ID(7, i, j, k)]*qsigsq);
          dwdsb[IDWS(2, 1, i, j, k)] = -qch * (qch*qch + sig*sig * (1 + prims[ID(6, i, j, k)]*prims[ID(6, i, j, k)]));
          dwdsb[IDWS(2, 2, i, j, k)] = -sig * (prims[ID(6, i, j, k)] * prims[ID(7, i, j, k)]*qch*sig - prims[ID(5, i, j, k)]*qsigsq);
          dwdsb[IDWS(3, 0, i, j, k)] = -sig * (prims[ID(5, i, j, k)] * prims[ID(7, i, j, k)]*qch*sig - prims[ID(6, i, j, k)]*qsigsq);
          dwdsb[IDWS(3, 1, i, j, k)] = -sig * (prims[ID(6, i, j, k)] * prims[ID(7, i, j, k)]*qch*sig + prims[ID(5, i, j, k)]*qsigsq);
          dwdsb[IDWS(3, 2, i, j, k)] = -qch * (qch*qch + sig*sig * (1 + prims[ID(7, i, j, k)]*prims[ID(7, i, j, k)]));
        }

        // Now do B
        {
          dwdsb[IDWS(5, 0, i, j, k)] = -prims[ID(5, i, j, k)] * sigcu * E[ID(0, i, j, k)];
          dwdsb[IDWS(5, 1, i, j, k)] = -prims[ID(6, i, j, k)] * sigcu * E[ID(0, i, j, k)] - sig*qsigsq*prims[ID(3, i, j, k)];
          dwdsb[IDWS(5, 2, i, j, k)] = -prims[ID(7, i, j, k)] * sigcu * E[ID(0, i, j, k)] + sig*qsigsq*prims[ID(2, i, j, k)];
          dwdsb[IDWS(6, 0, i, j, k)] = -prims[ID(5, i, j, k)] * sigcu * E[ID(1, i, j, k)] + sig*qsigsq*prims[ID(3, i, j, k)];
          dwdsb[IDWS(6, 1, i, j, k)] = -prims[ID(6, i, j, k)] * sigcu * E[ID(1, i, j, k)];
          dwdsb[IDWS(6, 2, i, j, k)] = -prims[ID(7, i, j, k)] * sigcu * E[ID(1, i, j, k)] - sig*qsigsq*prims[ID(1, i, j, k)];
          dwdsb[IDWS(7, 0, i, j, k)] = -prims[ID(5, i, j, k)] * sigcu * E[ID(2, i, j, k)] - sig*qsigsq*prims[ID(2, i, j, k)];
          dwdsb[IDWS(7, 1, i, j, k)] = -prims[ID(6, i, j, k)] * sigcu * E[ID(2, i, j, k)] + sig*qsigsq*prims[ID(1, i, j, k)];
          dwdsb[IDWS(7, 2, i, j, k)] = -prims[ID(7, i, j, k)] * sigcu * E[ID(2, i, j, k)];
        }

        // Now, C
        {
          dwdsb[IDWS(8, 0, i, j, k)] = -prims[ID(5, i, j, k)] * prims[ID(5, i, j, k)] * sigcu - sig*qsigsq;
          dwdsb[IDWS(8, 1, i, j, k)] = -prims[ID(5, i, j, k)] * prims[ID(6, i, j, k)] * sigcu;
          dwdsb[IDWS(8, 2, i, j, k)] = -prims[ID(5, i, j, k)] * prims[ID(7, i, j, k)] * sigcu;
          dwdsb[IDWS(9, 0, i, j, k)] = -prims[ID(6, i, j, k)] * prims[ID(5, i, j, k)] * sigcu;
          dwdsb[IDWS(9, 1, i, j, k)] = -prims[ID(6, i, j, k)] * prims[ID(6, i, j, k)] * sigcu - sig*qsigsq;
          dwdsb[IDWS(9, 2, i, j, k)] = -prims[ID(6, i, j, k)] * prims[ID(7, i, j, k)] * sigcu;
          dwdsb[IDWS(10, 0, i, j, k)] = -prims[ID(7, i, j, k)] * prims[ID(5, i, j, k)] * sigcu;
          dwdsb[IDWS(10, 1, i, j, k)] = -prims[ID(7, i, j, k)] * prims[ID(6, i, j, k)] * sigcu;
          dwdsb[IDWS(10, 2, i, j, k)] = -prims[ID(7, i, j, k)] * prims[ID(7, i, j, k)] * sigcu - sig*qsigsq;
        }

        // Finally, D
        {
          dwdsb[IDWS(11, 0, i, j, k)] = -prims[ID(1, i, j, k)]*qsigsq - sig*sig*prims[ID(5, i, j, k)]*vdotB;
          dwdsb[IDWS(11, 1, i, j, k)] = -prims[ID(2, i, j, k)]*qsigsq - sig*sig*prims[ID(6, i, j, k)]*vdotB;
          dwdsb[IDWS(11, 2, i, j, k)] = -prims[ID(3, i, j, k)]*qsigsq - sig*sig*prims[ID(7, i, j, k)]*vdotB;
        }

        // Dont forget to multiply by the prefactor! (if not testing)
        for (int l(0); l<12; l++) {
          for (int m(0); m<3; m++) {
            dwdsb[IDWS(l, m, i, j, k)] *= alpha[ID(0, i, j, k)];
          }
        }
      }
    }
  }
}

void ResistiveSGM::set_Dx(double * cons, double * prims, double * aux)
{
  // TODO
}

void ResistiveSGM::set_Dy(double * cons, double * prims, double * aux)
{
  // TODO
}


void ResistiveSGM::set_Dz(double * cons, double * prims, double * aux)
{
  // TODO
}


void ResistiveSGM::set_K(double * cons, double * prims, double * aux)
{
  // TODO
}


void ResistiveSGM::set_dfxdw(double * cons, double * prims, double * aux)
{
  // TODO
}


void ResistiveSGM::set_dfydw(double * cons, double * prims, double * aux)
{
  // TODO
}


void ResistiveSGM::set_dfzdw(double * cons, double * prims, double * aux)
{
  // TODO
}
