#include "gtest/gtest.h"
#include "platformEnv.h"
#include "boundaryConds.h"

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);

    PlatformEnv env = PlatformEnv(0, NULL, 1, 1, 1, 0);
    printf("env testing: %d\n", env.testing);
    return RUN_ALL_TESTS();
}
