#include "gtest/gtest.h"
#include "srmhd.h"
#include "simData.h"
#include "fluxVectorSplitting.h"
#include "resistiveSGM.h"

namespace
{

  TEST(RSGM, DataAssignment1D)
  /*
    Checks that, for 1-dimensional simulations, the variables are set correctly
  */
  {
    Data d(100, 0, 0, 0.01, 2.01, 0, 1, 0, 1, 0.4, 0.1, 4, 2, 50);
    SRMHD model(&d);
    FVS fluxMethod(&d, &model);
    ResistiveSGM subgridModel(&d, &fluxMethod);
  }


}
