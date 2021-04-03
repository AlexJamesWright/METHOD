// Serial main
#include "simData.h"
#include "simulation.h"
#include "initFunc.h"
#include "toy_q.h"
#include "toy_q_ce.h"
#include "boundaryConds.h"
#include "rkSplit.h"
#include "SSP2.h"
#include "fluxVectorSplitting.h"
#include "serialEnv.h"
#include "serialSaveDataHDF5.h"
#include "weno.h"
#include <cstring>

using namespace std;

int main(int argc, char *argv[]) {

  // Set up domain
  int Ng(4);
  int base_nx(16);
  int ny(0);
  int nz(0);
  double xmin(0.0);
  double xmax(1.0);
  double ymin(0.0);
  double ymax(1.0);
  double zmin(0.0);
  double zmax(1.0);
  double base_endTime(5.0);
  double cfl(0.4);
  bool output(false);
  int nreports(1);
  double base_kappa = 1.0e-3;
  double base_tau_q = 1.0e-3;

  SerialEnv env_imex(&argc, &argv, 1, 1, 1);
  SerialEnv env_ce(&argc, &argv, 1, 1, 1);

  std::string base_dir("1d/");
  std::string base_imex = base_dir+"data_imex";
  std::string base_ce = base_dir+"data_ce";

  for (int params_n=0; params_n<4; params_n++) {
    double kappa = base_kappa / pow(10, params_n);
    double tau_q = base_tau_q / pow(10, params_n);
    double endTime = base_endTime * pow(10, params_n);
    cout << "Doing params_n " << params_n << " " << kappa << endl;
    int nx = base_nx;
    while ((nx * kappa < 1.0) && (nx < 5000)) {
      cout << "Resolution " << nx << endl;

      DataArgs data_args_imex(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
      data_args_imex.sCfl(cfl);
      data_args_imex.sNg(Ng);
      const std::vector<double> toy_params { {kappa, tau_q} };
      const std::vector<std::string> toy_param_names = {"kappa", "tau_q"};
      const int n_toy_params(2);
      data_args_imex.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);

      Data data_imex(data_args_imex, &env_imex);

      DataArgs data_args_ce(nx, ny, nz, xmin, xmax, ymin, ymax, zmin, zmax, endTime);
      data_args_ce.sCfl(cfl);
      data_args_ce.sNg(Ng);
      // const std::vector<double> toy_params { {kappa, tau_q} };
      // const std::vector<std::string> toy_param_names = {"kappa", "tau_q"};
      // const int n_toy_params(2);
      data_args_ce.sOptionalSimArgs(toy_params, toy_param_names, n_toy_params);

      Data data_ce(data_args_ce, &env_ce);

      // Choose particulars of simulation
      ToyQ model_imex(&data_imex);
      ToyQ_CE model_ce(&data_ce);

      Weno3 weno_imex(&data_imex);
      Weno3 weno_ce(&data_ce);

      FVS fluxMethod_imex(&data_imex, &weno_imex, &model_imex);
      FVS fluxMethod_ce(&data_ce, &weno_ce, &model_ce);

      // Outflow bcs(&data);
      Periodic bcs_imex(&data_imex);
      Periodic bcs_ce(&data_ce);

      Simulation sim_imex(&data_imex, &env_imex);
      Simulation sim_ce(&data_ce, &env_ce);

      BlobToyQ init_imex(&data_imex);
      BlobToyQ_CE init_ce(&data_ce);
      // Blob2dToyQ init(&data);

      // RKSplit timeInt(&data, &model, &bcs, &fluxMethod);
      // BackwardsRK2 timeInt(&data, &model, &bcs, &fluxMethod);
      SSP2 timeInt_imex(&data_imex, &model_imex, &bcs_imex, &fluxMethod_imex);
      RKSplit timeInt_ce(&data_ce, &model_ce, &bcs_ce, &fluxMethod_ce);

      std::string params_string = "_nx"+std::to_string(nx)+"_kappa1em"+std::to_string(3+params_n);
      SerialSaveDataHDF5 save_imex(&data_imex, &env_imex, base_imex+params_string, SerialSaveDataHDF5::OUTPUT_ALL);
      SerialSaveDataHDF5 save_ce(&data_ce, &env_ce, base_ce+params_string, SerialSaveDataHDF5::OUTPUT_ALL);

      // Now objects have been created, set up the simulation
      sim_imex.set(&init_imex, &model_imex, &timeInt_imex, &bcs_imex, &fluxMethod_imex, &save_imex);
      sim_ce.set(&init_ce, &model_ce, &timeInt_ce, &bcs_ce, &fluxMethod_ce, &save_ce);

      sim_imex.evolve(output);
      save_imex.saveAll();
      sim_ce.evolve(output);
      save_ce.saveAll();

      nx *= 2;

    }
  }

  return 0;

}
