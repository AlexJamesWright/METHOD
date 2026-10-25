#include "subfilterGradient.h"
#include <cmath>
#include <stdexcept>

// Cell index within one variable; variables are stored as var * Nx*Ny*Nz + cell, as in ID(var, i, j, k)
#define CELL(i, j, k) (((i) * (d->Ny) + (j)) * (d->Nz) + (k))

SubfilterGradient::SubfilterGradient(Data * data, double c, int densityTerm, int limiter, double quadFactor) :
  ModelExtension(data), c(c), quadFactor(quadFactor), densityTerm(densityTerm), limiter(limiter), lastExchange(0.0), lastExchangeAbs(0.0)
{
  Data * d(this->data);
  if (d->dims != 3) throw std::invalid_argument("SubfilterGradient needs a three-dimensional simulation.\n");
  if (d->Ng < 2) throw std::invalid_argument("SubfilterGradient needs at least 2 ghost cells (gradient, then divergence).\n");
  if (d->Nprims < 8 || d->Ncons < 9) throw std::invalid_argument("SubfilterGradient is written for SRMHD (8 primitive, 9 conserved variables).\n");
  if (fabs(d->dx - d->dy) > 1e-12 * d->dx || fabs(d->dx - d->dz) > 1e-12 * d->dx) throw std::invalid_argument("SubfilterGradient needs cubic cells.\n");

  sourceExists = true;
  sigma2 = c * d->dx * d->dx;

  const long N(static_cast<long>(d->Nx) * d->Ny * d->Nz);
  E = new double[3 * N] ();
  Dens = new double[N] ();
  dv = new double[9 * N] ();
  dB = new double[9 * N] ();
  dE = new double[9 * N] ();
  dD = new double[3 * N] ();
  emf = new double[3 * N] ();
  poynting = new double[3 * N] ();
  stress = new double[9 * N] ();
  d->sourceExtension = new double[d->Ncons * N] ();
}

SubfilterGradient::~SubfilterGradient()
{
  Data * d(this->data);
  delete[] E;
  delete[] Dens;
  delete[] dv;
  delete[] dB;
  delete[] dE;
  delete[] dD;
  delete[] emf;
  delete[] poynting;
  delete[] stress;
  delete[] d->sourceExtension;
  d->sourceExtension = NULL;
}

void SubfilterGradient::sourceExtension(double * cons, double * prims, double * aux, double * source)
{
  Data * d(this->data);
  const int Nx(d->Nx), Ny(d->Ny), Nz(d->Nz);
  const long N(static_cast<long>(Nx) * Ny * Nz);
  const double sigma2q(quadFactor * sigma2);   // sigma^2 of the quadratic gradient parts of the stress and Poynting deficits
  const double inv2h[3] = {0.5 / d->dx, 0.5 / d->dy, 0.5 / d->dz};
  const long stride[3] = {static_cast<long>(Ny) * Nz, static_cast<long>(Nz), 1};

  for (long n(0); n < d->Ncons * N; n++) source[n] = 0.0;

  // Resolved electric field E_r = -v x B, and D = rho W, at every cell
  for (int i(0); i < Nx; i++) {
    for (int j(0); j < Ny; j++) {
      for (int k(0); k < Nz; k++) {
        const long cell(CELL(i, j, k));
        const double vx(prims[1 * N + cell]), vy(prims[2 * N + cell]), vz(prims[3 * N + cell]);
        const double Bx(prims[5 * N + cell]), By(prims[6 * N + cell]), Bz(prims[7 * N + cell]);
        E[0 * N + cell] = -(vy * Bz - vz * By);
        E[1 * N + cell] = -(vz * Bx - vx * Bz);
        E[2 * N + cell] = -(vx * By - vy * Bx);
        Dens[cell] = prims[0 * N + cell] / sqrt(1.0 - (vx * vx + vy * vy + vz * vz));
      }
    }
  }

  // Gradients, EMF, Poynting-flux and stress deficits at every cell with a full one-cell stencil
  for (int i(1); i < Nx - 1; i++) {
    for (int j(1); j < Ny - 1; j++) {
      for (int k(1); k < Nz - 1; k++) {
        const long cell(CELL(i, j, k));
        for (int l(0); l < 3; l++) {
          for (int m(0); m < 3; m++) {
            dv[(l * 3 + m) * N + cell] = (prims[(1 + m) * N + cell + stride[l]] - prims[(1 + m) * N + cell - stride[l]]) * inv2h[l];
            dB[(l * 3 + m) * N + cell] = (prims[(5 + m) * N + cell + stride[l]] - prims[(5 + m) * N + cell - stride[l]]) * inv2h[l];
          }
          dD[l * N + cell] = (Dens[cell + stride[l]] - Dens[cell - stride[l]]) * inv2h[l];
        }
        const double v[3] = {prims[1 * N + cell], prims[2 * N + cell], prims[3 * N + cell]};
        const double B[3] = {prims[5 * N + cell], prims[6 * N + cell], prims[7 * N + cell]};
        const double Er[3] = {E[0 * N + cell], E[1 * N + cell], E[2 * N + cell]};

        double dvl[3][3], dBl[3][3], dEl[3][3];
        for (int l(0); l < 3; l++) {
          for (int m(0); m < 3; m++) {
            dvl[l][m] = dv[(l * 3 + m) * N + cell];
            dBl[l][m] = dB[(l * 3 + m) * N + cell];
          }
          // d_l E = -d_l v x B - v x d_l B
          dEl[l][0] = -(dvl[l][1] * B[2] - dvl[l][2] * B[1]) - (v[1] * dBl[l][2] - v[2] * dBl[l][1]);
          dEl[l][1] = -(dvl[l][2] * B[0] - dvl[l][0] * B[2]) - (v[2] * dBl[l][0] - v[0] * dBl[l][2]);
          dEl[l][2] = -(dvl[l][0] * B[1] - dvl[l][1] * B[0]) - (v[0] * dBl[l][1] - v[1] * dBl[l][0]);
          for (int m(0); m < 3; m++) dE[(l * 3 + m) * N + cell] = dEl[l][m];
        }

        // EMF = -sigma^2 sum_l d_l v x d_l B  (+ zeta sigma^2 (sum_l d_l D d_l v / D) x B)
        double em[3] = {0.0, 0.0, 0.0};
        for (int l(0); l < 3; l++) {
          em[0] -= sigma2 * (dvl[l][1] * dBl[l][2] - dvl[l][2] * dBl[l][1]);
          em[1] -= sigma2 * (dvl[l][2] * dBl[l][0] - dvl[l][0] * dBl[l][2]);
          em[2] -= sigma2 * (dvl[l][0] * dBl[l][1] - dvl[l][1] * dBl[l][0]);
        }
        if (densityTerm) {
          double q[3] = {0.0, 0.0, 0.0};
          for (int l(0); l < 3; l++) for (int m(0); m < 3; m++) q[m] += dD[l * N + cell] * dvl[l][m] / Dens[cell];
          em[0] += sigma2 * (q[1] * B[2] - q[2] * B[1]);
          em[1] += sigma2 * (q[2] * B[0] - q[0] * B[2]);
          em[2] += sigma2 * (q[0] * B[1] - q[1] * B[0]);
        }
        for (int m(0); m < 3; m++) emf[m * N + cell] = em[m];

        // Poynting-flux deficit = EMF x B + sigma^2 sum_l d_l E x d_l B
        double po[3] = {em[1] * B[2] - em[2] * B[1], em[2] * B[0] - em[0] * B[2], em[0] * B[1] - em[1] * B[0]};
        for (int l(0); l < 3; l++) {
          po[0] += sigma2q * (dEl[l][1] * dBl[l][2] - dEl[l][2] * dBl[l][1]);
          po[1] += sigma2q * (dEl[l][2] * dBl[l][0] - dEl[l][0] * dBl[l][2]);
          po[2] += sigma2q * (dEl[l][0] * dBl[l][1] - dEl[l][1] * dBl[l][0]);
        }
        for (int m(0); m < 3; m++) poynting[m * N + cell] = po[m];

        // Maxwell-stress deficit
        const double Eem(Er[0] * em[0] + Er[1] * em[1] + Er[2] * em[2]);
        double quadTrace(0.0);
        for (int l(0); l < 3; l++) {
          for (int m(0); m < 3; m++) quadTrace += dEl[l][m] * dEl[l][m] + dBl[l][m] * dBl[l][m];
        }
        for (int a(0); a < 3; a++) {
          for (int b(0); b < 3; b++) {
            double s(-(Er[a] * em[b] + em[a] * Er[b]));
            if (a == b) s += Eem + 0.5 * sigma2q * quadTrace;
            for (int l(0); l < 3; l++) s += sigma2q * (-dEl[l][a] * dEl[l][b] - dBl[l][a] * dBl[l][b]);
            stress[(a * 3 + b) * N + cell] = s;
          }
        }
      }
    }
  }

  // Sources over the interior: -curl EMF (B), -div stress (S), -div Poynting (tau)
  lastExchange = 0.0;
  lastExchangeAbs = 0.0;
  const double vol(d->dx * d->dy * d->dz);
  for (int i(d->is); i < d->ie; i++) {
    for (int j(d->js); j < d->je; j++) {
      for (int k(d->ks); k < d->ke; k++) {
        const long cell(CELL(i, j, k));
        double dem[3][3];   // dem[l][m] = d_l EMF_m
        for (int l(0); l < 3; l++) {
          for (int m(0); m < 3; m++) dem[l][m] = (emf[m * N + cell + stride[l]] - emf[m * N + cell - stride[l]]) * inv2h[l];
        }
        const double sB[3] = {-(dem[1][2] - dem[2][1]), -(dem[2][0] - dem[0][2]), -(dem[0][1] - dem[1][0])};
        double sS[3], sTau(0.0);
        for (int a(0); a < 3; a++) {
          sS[a] = 0.0;
          for (int l(0); l < 3; l++) sS[a] -= (stress[(a * 3 + l) * N + cell + stride[l]] - stress[(a * 3 + l) * N + cell - stride[l]]) * inv2h[l];
        }
        for (int l(0); l < 3; l++) sTau -= (poynting[l * N + cell + stride[l]] - poynting[l * N + cell - stride[l]]) * inv2h[l];

        // Local energy exchange J.EMF with J = curl B
        const double J[3] = {dB[(1 * 3 + 2) * N + cell] - dB[(2 * 3 + 1) * N + cell],
                             dB[(2 * 3 + 0) * N + cell] - dB[(0 * 3 + 2) * N + cell],
                             dB[(0 * 3 + 1) * N + cell] - dB[(1 * 3 + 0) * N + cell]};
        const double exchange(J[0] * emf[0 * N + cell] + J[1] * emf[1 * N + cell] + J[2] * emf[2 * N + cell]);
        lastExchange += exchange * vol;
        lastExchangeAbs += fabs(exchange) * vol;
        if (limiter == 1 && exchange < 0.0) continue;

        for (int a(0); a < 3; a++) {
          source[(1 + a) * N + cell] = sS[a];
          source[(5 + a) * N + cell] = sB[a];
        }
        source[4 * N + cell] = sTau;
      }
    }
  }
}
