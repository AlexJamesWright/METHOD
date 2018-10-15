#include "resistiveSGM.h"
#include <cstdio>
#include <cmath>

#define ID(variable, idx, jdx, kdx)  ((variable)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

// dwdsb
#define IDWS(ldx, mdx, idx, jdx, kdx)  ((ldx)*(3)*(d->Nx)*(d->Ny)*(d->Nz) + (mdx)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
// dfxdw, dfydw, dfzdw
#define IDFW(ldx, mdx, idx, jdx, kdx)  ((ldx)*(12)*(d->Nx)*(d->Ny)*(d->Nz) + (mdx)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))
// Mx, My, and Mz matrix
#define IDM(ldx, mdx, idx, jdx, kdx)  ((ldx)*(3)*(d->Nx)*(d->Ny)*(d->Nz) + (mdx)*(d->Nx)*(d->Ny)*(d->Nz) + (idx)*(d->Ny)*(d->Nz) + (jdx)*(d->Nz) + (kdx))

void minmod(double * array, double * ahead, double * behind, double * source, double delta, int Nv, int Nx, int Ny, int Nz, int dir, Data * d);

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
  delete[] dfxdw;
  delete[] dfydw;
  delete[] dfzdw;
  delete[] dwdsb;
  delete[] E;
  delete[] q;
  delete[] K; //
  delete[] Mx; //
  delete[] My; //
  delete[] Mz; //
  delete[] fx;
  delete[] fy;
  delete[] fz;
  delete[] diffuX; //
  delete[] diffuY; //
  delete[] diffuZ; //
  delete[] alpha;
}

void ResistiveSGM::subgridSource(double * cons, double * prims, double * aux, double * source)
{
  // Syntax
  Data * d(this->data);

  // Zero arrays and set vars
  this->set_vars(cons, prims, aux);
  this->reset(source);

  // Ensure K and dwdsb are set
  this->set_K(cons, prims, aux);
  this->set_dwdsb(cons, prims, aux);

  // MINMOD
  // For this, I will use data->f and data->fnet as the ahead and behind
  // arrays, but will rename, saving on memory.
  double * ahead(d->f);
  double * behind(d->fnet);

  // Do minmod on Da to get gradient and add to source
  {
    // Determine the diffusion vectors
    this->set_Dx(cons, prims, aux);
    minmod(diffuX, ahead, behind, source, d->dx, d->Ncons, d->Nx, d->Ny, d->Nz, 0, d);
    if (d->dims>1)
    {
      this->set_Dy(cons, prims, aux);
      minmod(diffuY, ahead, behind, source, d->dy, d->Ncons, d->Nx, d->Ny, d->Nz, 1, d);
    }
    if (d->dims==3)
    {
      this->set_Dz(cons, prims, aux);
      minmod(diffuZ, ahead, behind, source, d->dz, d->Ncons, d->Nx, d->Ny, d->Nz, 2, d);
    }
  }
}

void ResistiveSGM::reset(double * source)
{
  // Syntax
  Data * d(this->data);
  // Reset the arrays in which we use the += operator
  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        // Source vector, and Da
        for (int var(0); var<d->Ncons; var++) {
          source[ID(var, i, j, k)] = 0.0;
        }
        for (int var(0); var<9; var++) {
          diffuX[ID(var, i, j, k)] = 0.0;
          diffuY[ID(var, i, j, k)] = 0.0;
          diffuZ[ID(var, i, j, k)] = 0.0;
        }
        // partial_a f^a
        K[ID(0, i, j, k)] = 0.0;
        K[ID(1, i, j, k)] = 0.0;
        K[ID(2, i, j, k)] = 0.0;
        // Mx, My, Mz
        for (int l(0); l<9; l++) {
          for (int m(0); m<3; m++) {
            Mx[IDM(l, m, i, j, k)] = 0.0;
            My[IDM(l, m, i, j, k)] = 0.0;
            Mz[IDM(l, m, i, j, k)] = 0.0;
          }
        }
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

        // Dont forget to multiply by the prefactor!
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
  // Syntax
  Data * d(this->data);

  this->set_dfxdw(cons, prims, aux);
  // Mx = -1 * DOT(dfxdw, dwdsb)
  for (int l(0); l<9; l++) {
    for (int m(0); m<3; m++) {
      for (int i(0); i<d->Nx; i++) {
        for (int j(0); j<d->Ny; j++) {
          for (int k(0); k<d->Nz; k++) {
            for (int n(0); n<12; n++) {
              Mx[IDM(l, m, i, j, k)] -= dfxdw[IDFW(l, n, i, j, k)] * dwdsb[IDWS(n, m, i, j, k)];
            }
          }
        }
      }
    }
  }
  // Dx = DOT(Mx, K)
  for (int l(0); l<9; l++) {
    for (int i(0); i<d->Nx; i++) {
      for (int j(0); j<d->Ny; j++) {
        for (int k(0); k<d->Nz; k++) {
          for (int m(0); m<3; m++) {
            diffuX[ID(l, i, j, k)] += Mx[IDM(l, m, i, j, k)] * K[ID(m, i, j, k)];
          }
        }
      }
    }
  }
}

void ResistiveSGM::set_Dy(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  this->set_dfydw(cons, prims, aux);


  // My = -1 * DOT(dfydw, dwdsb)
  for (int l(0); l<9; l++) {
    for (int m(0); m<3; m++) {
      for (int i(0); i<d->Nx; i++) {
        for (int j(0); j<d->Ny; j++) {
          for (int k(0); k<d->Nz; k++) {
            for (int n(0); n<12; n++) {
              My[IDM(l, m, i, j, k)] -= dfydw[IDFW(l, n, i, j, k)] * dwdsb[IDWS(n, m, i, j, k)];
            }
          }
        }
      }
    }
  }
  // Dy = DOT(My, K)
  for (int l(0); l<9; l++) {
    for (int i(0); i<d->Nx; i++) {
      for (int j(0); j<d->Ny; j++) {
        for (int k(0); k<d->Nz; k++) {
          for (int m(0); m<3; m++) {
            diffuY[ID(l, i, j, k)] += My[IDM(l, m, i, j, k)] * K[ID(m, i, j, k)];
          }
        }
      }
    }
  }
}


void ResistiveSGM::set_Dz(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  this->set_dfzdw(cons, prims, aux);
  // Mz = -1 * DOT(dfzdw, dwdsb)
  for (int l(0); l<9; l++) {
    for (int m(0); m<3; m++) {
      for (int i(0); i<d->Nx; i++) {
        for (int j(0); j<d->Ny; j++) {
          for (int k(0); k<d->Nz; k++) {
            for (int n(0); n<12; n++) {
              Mz[IDM(l, m, i, j, k)] -= dfzdw[IDFW(l, n, i, j, k)] * dwdsb[IDWS(n, m, i, j, k)];
            }
          }
        }
      }
    }
  }
  // Dz = DOT(Mz, K)
  for (int l(0); l<9; l++) {
    for (int i(0); i<d->Nx; i++) {
      for (int j(0); j<d->Ny; j++) {
        for (int k(0); k<d->Nz; k++) {
          for (int m(0); m<3; m++) {
            diffuZ[ID(l, i, j, k)] += Mz[IDM(l, m, i, j, k)] * K[ID(m, i, j, k)];
          }
        }
      }
    }
  }
}


void ResistiveSGM::set_K(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

// First do x-direction
  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        fx[ID(0, i, j, k)] = 0;
        fx[ID(1, i, j, k)] = prims[ID(7, i, j, k)];
        fx[ID(2, i, j, k)] = -prims[ID(6, i, j, k)];
      }
    }
  }
  // Reconstruct stiff fluxes
  this->fluxMethod->fluxReconstruction(E, NULL, NULL, fx, d->fnet, 0, 3);
  // Add flux differencing to K
  for (int var(0); var<3; var++) {
    for (int i(1); i<d->Nx-1; i++) {
      for (int j(0); j<d->Ny; j++) {
        for (int k(0); k<d->Nz; k++) {
          K[ID(var, i, j, k)] += d->fnet[ID(var, i+1, j, k)]/(2*d->dx) - d->fnet[ID(var, i-1, j, k)]/(2*d->dx);
        }
      }
    }
  }

  // Now add y-contribution
  if (d->dims>1)
  {
    for (int i(0); i<d->Nx; i++) {
      for (int j(0); j<d->Ny; j++) {
        for (int k(0); k<d->Nz; k++) {
          fy[ID(0, i, j, k)] = -prims[ID(7, i, j, k)];
          fy[ID(1, i, j, k)] = 0;
          fy[ID(2, i, j, k)] = prims[ID(5, i, j, k)];
        }
      }
    }
    // Reconstruct stiff fluxes
    this->fluxMethod->fluxReconstruction(E, NULL, NULL, fy, d->fnet, 1, 3);
    // Add flux differencing to K
    for (int var(0); var<3; var++) {
      for (int i(0); i<d->Nx; i++) {
        for (int j(1); j<d->Ny-1; j++) {
          for (int k(0); k<d->Nz; k++) {
            K[ID(var, i, j, k)] += d->fnet[ID(var, i, j+1, k)]/(2*d->dy) - d->fnet[ID(var, i, j-1, k)]/(2*d->dy);
          }
        }
      }
    }
  }

  // Finally, add z-contribution
  if (d->dims==3)
  {
    for (int i(0); i<d->Nx; i++) {
      for (int j(0); j<d->Ny; j++) {
        for (int k(0); k<d->Nz; k++) {
          fz[ID(0, i, j, k)] = prims[ID(6, i, j, k)];
          fz[ID(1, i, j, k)] = -prims[ID(5, i, j, k)];
          fz[ID(2, i, j, k)] = 0;
        }
      }
    }
    // Reconstruct stiff fluxes
    this->fluxMethod->fluxReconstruction(E, NULL, NULL, fz, d->fnet, 2, 3);
    // Add flux differencing to K
    for (int var(0); var<3; var++) {
      for (int i(0); i<d->Nx; i++) {
        for (int j(0); j<d->Ny; j++) {
          for (int k(1); k<d->Nz-1; k++) {
            K[ID(var, i, j, k)] += d->fnet[ID(var, i, j, k+1)]/(2*d->dz) - d->fnet[ID(var, i, j, k-1)]/(2*d->dz);
          }
        }
      }
    }

  }

}


void ResistiveSGM::set_dfxdw(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        // Row 0
        dfxdw[IDFW(0, 0, i, j, k)] = prims[ID(1, i, j, k)];
        dfxdw[IDFW(0, 1, i, j, k)] = prims[ID(0, i, j, k)];
        // Row 1
        dfxdw[IDFW(1, 4, i, j, k)] = 1;
        dfxdw[IDFW(1, 5, i, j, k)] = -prims[ID(5, i, j, k)];
        dfxdw[IDFW(1, 6, i, j, k)] = prims[ID(6, i, j, k)];
        dfxdw[IDFW(1, 7, i, j, k)] = prims[ID(7, i, j, k)];
        dfxdw[IDFW(1, 8, i, j, k)] = -E[ID(0, i, j, k)];
        dfxdw[IDFW(1, 9, i, j, k)] = E[ID(1, i, j, k)];
        dfxdw[IDFW(1, 10, i, j, k)] = E[ID(2, i, j, k)];
        // Row 2
        dfxdw[IDFW(2, 5, i, j, k)] = -prims[ID(6, i, j, k)];
        dfxdw[IDFW(2, 6, i, j, k)] = -prims[ID(5, i, j, k)];
        dfxdw[IDFW(2, 8, i, j, k)] = -E[ID(1, i, j, k)];
        dfxdw[IDFW(2, 9, i, j, k)] = -E[ID(0, i, j, k)];
        // Row 3
        dfxdw[IDFW(3, 5, i, j, k)] = -prims[ID(7, i, j, k)];
        dfxdw[IDFW(3, 7, i, j, k)] = -prims[ID(5, i, j, k)];
        dfxdw[IDFW(3, 8, i, j, k)] = -E[ID(2, i, j, k)];
        dfxdw[IDFW(3, 10, i, j, k)] = -E[ID(0, i, j, k)];
        // Row 4
        dfxdw[IDFW(4, 1, i, j, k)] = d->gamma*prims[ID(4, i, j, k)] / (d->gamma - 1);
        dfxdw[IDFW(4, 4, i, j, k)] = d->gamma*prims[ID(1, i, j, k)] / (d->gamma - 1);
        dfxdw[IDFW(4, 6, i, j, k)] = -E[ID(2, i, j, k)];
        dfxdw[IDFW(4, 7, i, j, k)] = E[ID(1, i, j, k)];
        dfxdw[IDFW(4, 9, i, j, k)] = prims[ID(7, i, j, k)];
        dfxdw[IDFW(4, 10, i, j, k)] = -prims[ID(6, i, j, k)];
        // Row 6
        dfxdw[IDFW(6, 10, i, j, k)] = -1;
        // Row 7
        dfxdw[IDFW(7, 9, i, j, k)] = 1;
        // Row 8
        dfxdw[IDFW(8, 1, i, j, k)] = q[ID(0, i, j, k)];
        dfxdw[IDFW(8, 2, i, j, k)] = d->sigma*prims[ID(7, i, j, k)];
        dfxdw[IDFW(8, 3, i, j, k)] = -d->sigma*prims[ID(6, i, j, k)];
        dfxdw[IDFW(8, 6, i, j, k)] = -d->sigma*prims[ID(3, i, j, k)];
        dfxdw[IDFW(8, 7, i, j, k)] = d->sigma*prims[ID(2, i, j, k)];
        dfxdw[IDFW(8, 8, i, j, k)] = d->sigma;
        dfxdw[IDFW(8, 11, i, j, k)] = prims[ID(1, i, j, k)];

      }
    }
  }
}


void ResistiveSGM::set_dfydw(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        // Row 0
        dfydw[IDFW(0, 0, i, j, k)] = prims[ID(2, i, j, k)];
        dfydw[IDFW(0, 2, i, j, k)] = prims[ID(0, i, j, k)];
        // Row 1
        dfydw[IDFW(1, 5, i, j, k)] = -prims[ID(6, i, j, k)];
        dfydw[IDFW(1, 6, i, j, k)] = -prims[ID(5, i, j, k)];
        dfydw[IDFW(1, 8, i, j, k)] = -E[ID(1, i, j, k)];
        dfydw[IDFW(1, 9, i, j, k)] = -E[ID(0, i, j, k)];
        // Row 2
        dfydw[IDFW(2, 4, i, j, k)] = 1;
        dfydw[IDFW(2, 5, i, j, k)] = prims[ID(5, i, j, k)];
        dfydw[IDFW(2, 6, i, j, k)] = -prims[ID(6, i, j, k)];
        dfydw[IDFW(2, 7, i, j, k)] = prims[ID(7, i, j, k)];
        dfydw[IDFW(2, 8, i, j, k)] = E[ID(0, i, j, k)];
        dfydw[IDFW(2, 9, i, j, k)] = -E[ID(1, i, j, k)];
        dfydw[IDFW(2, 10, i, j, k)] = E[ID(2, i, j, k)];
        // Row 3
        dfydw[IDFW(3, 6, i, j, k)] = -prims[ID(7, i, j, k)];
        dfydw[IDFW(3, 7, i, j, k)] = -prims[ID(6, i, j, k)];
        dfydw[IDFW(3, 9, i, j, k)] = -E[ID(2, i, j, k)];
        dfydw[IDFW(3, 10, i, j, k)] = -E[ID(1, i, j, k)];
        // Row 4
        dfydw[IDFW(4, 2, i, j, k)] = d->gamma*prims[ID(4, i, j, k)] / (d->gamma - 1);
        dfydw[IDFW(4, 4, i, j, k)] = d->gamma*prims[ID(2, i, j, k)] / (d->gamma - 1);
        dfydw[IDFW(4, 5, i, j, k)] = E[ID(2, i, j, k)];
        dfydw[IDFW(4, 7, i, j, k)] = -E[ID(0, i, j, k)];
        dfydw[IDFW(4, 8, i, j, k)] = -prims[ID(7, i, j, k)];
        dfydw[IDFW(4, 10, i, j, k)] = prims[ID(5, i, j, k)];
        // Row 6
        dfydw[IDFW(5, 10, i, j, k)] = 1;
        // Row 7
        dfydw[IDFW(7, 8, i, j, k)] = -1;
        // Row 8
        dfydw[IDFW(8, 1, i, j, k)] = -d->sigma*prims[ID(7, i, j, k)];
        dfydw[IDFW(8, 2, i, j, k)] = q[ID(0, i, j, k)];
        dfydw[IDFW(8, 3, i, j, k)] = d->sigma*prims[ID(5, i, j, k)];
        dfydw[IDFW(8, 5, i, j, k)] = d->sigma*prims[ID(3, i, j, k)];
        dfydw[IDFW(8, 7, i, j, k)] = -d->sigma*prims[ID(1, i, j, k)];
        dfydw[IDFW(8, 9, i, j, k)] = d->sigma;
        dfydw[IDFW(8, 11, i, j, k)] = prims[ID(2, i, j, k)];

      }
    }
  }
}


void ResistiveSGM::set_dfzdw(double * cons, double * prims, double * aux)
{
  // Syntax
  Data * d(this->data);

  for (int i(0); i<d->Nx; i++) {
    for (int j(0); j<d->Ny; j++) {
      for (int k(0); k<d->Nz; k++) {
        // Row 0
        dfzdw[IDFW(0, 0, i, j, k)] = prims[ID(3, i, j, k)];
        dfzdw[IDFW(0, 3, i, j, k)] = prims[ID(0, i, j, k)];
        // Row 1
        dfzdw[IDFW(1, 5, i, j, k)] = -prims[ID(7, i, j, k)];
        dfzdw[IDFW(1, 7, i, j, k)] = -prims[ID(5, i, j, k)];
        dfzdw[IDFW(1, 8, i, j, k)] = -E[ID(2, i, j, k)];
        dfzdw[IDFW(1, 10, i, j, k)] = -E[ID(0, i, j, k)];
        // Row 2
        dfzdw[IDFW(2, 6, i, j, k)] = -prims[ID(7, i, j, k)];
        dfzdw[IDFW(2, 7, i, j, k)] = -prims[ID(6, i, j, k)];
        dfzdw[IDFW(2, 9, i, j, k)] = -E[ID(2, i, j, k)];
        dfzdw[IDFW(2, 10, i, j, k)] = -E[ID(1, i, j, k)];
        // Row 3
        dfzdw[IDFW(3, 4, i, j, k)] = 1;
        dfzdw[IDFW(3, 5, i, j, k)] = prims[ID(5, i, j, k)];
        dfzdw[IDFW(3, 6, i, j, k)] = prims[ID(6, i, j, k)];
        dfzdw[IDFW(3, 7, i, j, k)] = -prims[ID(7, i, j, k)];
        dfzdw[IDFW(3, 8, i, j, k)] = E[ID(0, i, j, k)];
        dfzdw[IDFW(3, 9, i, j, k)] = E[ID(1, i, j, k)];
        dfzdw[IDFW(3, 10, i, j, k)] = -E[ID(2, i, j, k)];
        // Row 4
        dfzdw[IDFW(4, 3, i, j, k)] = d->gamma*prims[ID(4, i, j, k)] / (d->gamma - 1);
        dfzdw[IDFW(4, 4, i, j, k)] = d->gamma*prims[ID(3, i, j, k)] / (d->gamma - 1);
        dfzdw[IDFW(4, 5, i, j, k)] = -E[ID(1, i, j, k)];
        dfzdw[IDFW(4, 6, i, j, k)] = E[ID(0, i, j, k)];
        dfzdw[IDFW(4, 8, i, j, k)] = prims[ID(6, i, j, k)];
        dfzdw[IDFW(4, 9, i, j, k)] = -prims[ID(5, i, j, k)];
        // Row 6
        dfzdw[IDFW(5, 9, i, j, k)] = -1;
        // Row 7
        dfzdw[IDFW(6, 8, i, j, k)] = 1;
        // Row 8
        dfzdw[IDFW(8, 1, i, j, k)] = d->sigma*prims[ID(6, i, j, k)];
        dfzdw[IDFW(8, 2, i, j, k)] = - d->sigma*prims[ID(5, i, j, k)];
        dfzdw[IDFW(8, 3, i, j, k)] = q[ID(0, i, j, k)];
        dfzdw[IDFW(8, 5, i, j, k)] = -d->sigma*prims[ID(2, i, j, k)];
        dfzdw[IDFW(8, 6, i, j, k)] = d->sigma*prims[ID(1, i, j, k)];
        dfzdw[IDFW(8, 10, i, j, k)] = d->sigma;
        dfzdw[IDFW(8, 11, i, j, k)] = prims[ID(3, i, j, k)];

      }
    }
  }
}


void minmod(double * array, double * ahead, double * behind, double * source, double delta, int Nv, int Nx, int Ny, int Nz, int dir, Data * d)
{

  /* ##################################################################
      This can be optimised.
      I do not need to store all ahead and behind. Store in a single variable at a time
      and calculate as needed.
      ################################################################*/

  if (dir == 0)
  {
    for (int var(0); var<Nv; var++) {
      for (int i(0); i<d->Nx-1; i++) {
        for (int j(0); j<d->Ny; j++) {
          for (int k(0);k<d->Nz; k++) {
            ahead[ID(var, i, j, k)] = (array[ID(var, i+1, j, k)] - array[ID(var, i, j, k)]) / delta;
            behind[ID(var, i+1, j, k)] = ahead[ID(var, i, j, k)];
          }
        }
      }
      for (int i(0); i<d->Nx-1; i++) {
        for (int j(0); j<d->Ny; j++) {
          for (int k(0);k<d->Nz; k++) {
            if (ahead[ID(var, i, j, k)]*behind[ID(var, i, j, k)] > 0) {
              if (fabs(ahead[ID(var, i, j, k)]) < fabs(behind[ID(var, i, j, k)]))
                source[ID(var, i, j, k)] += ahead[ID(var, i, j, k)];
              else
                source[ID(var, i, j, k)] += behind[ID(var, i, j, k)];
            }
          }
        }
      }
    }
  }


  else if (dir == 1)
  {
    for (int var(0); var<Nv; var++) {
      for (int i(0); i<d->Nx; i++) {
        for (int j(0); j<d->Ny-1; j++) {
          for (int k(0);k<d->Nz; k++) {
            ahead[ID(var, i, j, k)] = (array[ID(var, i, j+1, k)] - array[ID(var, i, j, k)]) / delta;
            behind[ID(var, i, j+1, k)] = ahead[ID(var, i, j, k)];
          }
        }
      }
      for (int i(0); i<d->Nx; i++) {
        for (int j(0); j<d->Ny-1; j++) {
          for (int k(0);k<d->Nz; k++) {
            if (ahead[ID(var, i, j, k)]*behind[ID(var, i, j, k)] > 0) {
              if (fabs(ahead[ID(var, i, j, k)]) < fabs(behind[ID(var, i, j, k)]))
                source[ID(var, i, j, k)] += ahead[ID(var, i, j, k)];
              else
                source[ID(var, i, j, k)] += behind[ID(var, i, j, k)];
            }
          }
        }
      }
    }
  }


  else
  {
    for (int var(0); var<Nv; var++) {
      for (int i(0); i<d->Nx; i++) {
        for (int j(0); j<d->Ny; j++) {
          for (int k(0); k<d->Nz-1; k++) {
            ahead[ID(var, i, j, k)] = (array[ID(var, i, j, k+1)] - array[ID(var, i, j, k)]) / delta;
            behind[ID(var, i, j, k+1)] = ahead[ID(var, i, j, k)];
          }
        }
      }
      for (int i(0); i<d->Nx; i++) {
        for (int j(0); j<d->Ny; j++) {
          for (int k(0);k<d->Nz-1; k++) {
            if (ahead[ID(var, i, j, k)]*behind[ID(var, i, j, k)] > 0) {
              if (fabs(ahead[ID(var, i, j, k)]) < fabs(behind[ID(var, i, j, k)]))
                source[ID(var, i, j, k)] += ahead[ID(var, i, j, k)];
              else
                source[ID(var, i, j, k)] += behind[ID(var, i, j, k)];
            }
          }
        }
      }
    }
  }
}
