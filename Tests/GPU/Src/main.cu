#include "gtest/gtest.h"
#include "parallelEnv.h"

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    // Create env here to ensure MPI initialisation is handled. Will need to create this object again inside each test
    // -- mpi init will only be called the first time
    ParallelEnv env(0, NULL, 2, 2, 1);
    return RUN_ALL_TESTS();
}
