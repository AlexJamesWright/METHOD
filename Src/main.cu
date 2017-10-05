#include "simulation.h"
#include "model.h"

int main(void) {

  Simulation sim(10, 10, 0, 1, 0, 1, 1);
  SRMHD model(&sim);

  return 0;

}
