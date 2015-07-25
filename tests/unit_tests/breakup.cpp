#include <gtest/gtest.h>
#include "breakupKernels/binaryBreakup.H"

TEST(breakupKernel, source){
    Foam::dimensionedScalar xi(
        "xi", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::dimensionedScalar expected(
        "expected", Foam::dimensionSet(0, 0, -1, 0, 0), 36);

    Foam::breakupKernels::binaryBreakup kernel;

    auto result = kernel.S(xi);

    ASSERT_EQ(expected.value(), result.value());
    ASSERT_EQ(expected.dimensions(), result.dimensions());
}
