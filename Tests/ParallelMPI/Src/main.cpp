#include "gtest/gtest.h"
#include "platformEnv.h"
#include "boundaryConds.h"

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    // Create env here to ensure MPI initialisation is handled. Will need to create this object again inside each test
    // -- in other calls the testing flag to platformEnv will be used to avoid calling mpi init/finalize
    PlatformEnv env = PlatformEnv(0, NULL, 1, 1, 1, 0);
    printf("env testing: %d\n", env.testing);
    return RUN_ALL_TESTS();
}
