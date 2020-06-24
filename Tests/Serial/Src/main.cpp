#include "gtest/gtest.h"
#include "platformEnv.h"

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    // Create env here to ensure MPI initialisation is handled. Will need to create this object again inside each test
    // -- mpi init will only be called the first time
    PlatformEnv env = PlatformEnv(0, NULL, 1, 1, 1);
    return RUN_ALL_TESTS();
}
