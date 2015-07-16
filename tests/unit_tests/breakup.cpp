#include <gtest/gtest.h>
#include "breakupKernels/binaryBreakup.H"

TEST(GSL_integrator, source){
    Foam::dimensionedScalar xi(
        "a_class", Foam::dimensionSet(0, 3, 0, 0, 0), 6);
    Foam::dimensionedScalar result(
        "result", Foam::dimensionSet(0, 3, 0, 0, 0), 36);
    Foam::breakupKernels::binaryBreakup kernel;
    ASSERT_EQ(result.value(), kernel.S(xi).value());
}
