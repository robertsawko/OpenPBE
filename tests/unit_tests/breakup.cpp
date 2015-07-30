#include <gtest/gtest.h>
#include "breakupKernels/binaryBreakup.H"
#include "breakupKernels/noBreakup.H"
#include "breakupKernels/CoulaloglouTavlarides.H"

TEST(breakup, binary){
    Foam::dimensionedScalar xi(
        "xi", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::dimensionedScalar expected(
        "expected", Foam::dimensionSet(0, 0, -1, 0, 0), 36);

    Foam::breakupKernels::binaryBreakupImpl kernel;

    auto result = kernel.S(xi);

    ASSERT_EQ(expected.value(), result.value());
    ASSERT_EQ(expected.dimensions(), result.dimensions());
}

TEST(breakup, none){
    Foam::dimensionedScalar xi(
        "xi", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::dimensionedScalar expected(
        "expected", Foam::dimensionSet(0, 0, -1, 0, 0), 0);

    Foam::breakupKernels::noBreakupImpl kernel;

    auto result = kernel.S(xi);

    ASSERT_EQ(expected.value(), result.value());
    ASSERT_EQ(expected.dimensions(), result.dimensions());
}

TEST(breakup, CoulaloglouTavlarides){
    Foam::scalar C1(0.0152), C2(0.0678),
                 sigma(0.04282), alpha(0.050000000000000003);

    Foam::dimensionedScalar epsilon("epsilon", Foam::dimless,
                                    0.12935361098784559);
    Foam::dimensionedScalar rho_d("rho_d", Foam::dimless, 972.0);

    Foam::dimensionedScalar xi(
        "xi", Foam::dimensionSet(0, 3, 0, 0, 0), 9.738310380055584507e-11);

    Foam::dimensionedScalar expected(
        "expected", Foam::dimensionSet(0, 0, -1, 0, 0),
                    1.122262582852271567e-02);

    Foam::breakupKernels::CoulaloglouTavlaridesImp kernel(C1, C2, alpha, sigma);

    auto result = kernel.S(xi, rho_d, epsilon);

    ASSERT_EQ(expected.value(), result.value());
    ASSERT_EQ(expected.dimensions(), result.dimensions());
}
