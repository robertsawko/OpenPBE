#include <gtest/gtest.h>
#include "coalescenceKernels/constantCoalescence.H"
#include "coalescenceKernels/noCoalescence.H"
#include "coalescenceKernels/CoulaloglouTavlarides.H"

TEST(coalescence, constant){
    Foam::dimensionedScalar xi1(
        "xi1", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::dimensionedScalar xi2(
        "xi2", Foam::dimensionSet(0, 3, 0, 0, 0), 8);

    Foam::scalar c = 5;

    Foam::coalescenceKernels::constantCoalescenceImpl kernel(c);

    Foam::dimensionedScalar expected(
        "expected", Foam::dimensionSet(0, 0, -1, 0, 0), c);

    auto result = kernel.S(xi1, xi2);

    ASSERT_EQ(expected.value(), result.value());
    ASSERT_EQ(expected.dimensions(), result.dimensions());
}

TEST(coalescence, none){
    Foam::dimensionedScalar xi1(
        "xi1", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::dimensionedScalar xi2(
        "xi2", Foam::dimensionSet(0, 3, 0, 0, 0), 8);

    Foam::coalescenceKernels::noCoalescenceImpl kernel;

    Foam::dimensionedScalar expected(
        "expected", Foam::dimensionSet(0, 0, -1, 0, 0), 0);

    auto result = kernel.S(xi1, xi2);

    ASSERT_EQ(expected.value(), result.value());
    ASSERT_EQ(expected.dimensions(), result.dimensions());
}

TEST(coalescence, CoulaloglouTavlarides){
    Foam::scalar C1(0.0152), C2(0.0678),
                 sigma(0.04282), alpha(0.050000000000000003);

    Foam::dimensionedScalar epsilon("epsilon", Foam::dimless,
                                    0.12935361098784559);
    Foam::dimensionedScalar rho_d("rho_d", Foam::dimless, 972.0);

    Foam::dimensionedScalar nud("nu_d", Foam::dimless, 0.001); //actually muc

    Foam::dimensionedScalar xi1(
        "xi1", Foam::dimensionSet(0, 3, 0, 0, 0), 5.071423503403169858e-10);

    Foam::dimensionedScalar xi2(
        "xi2", Foam::dimensionSet(0, 3, 0, 0, 0), 2.885425297794245861e-11);

    Foam::dimensionedScalar expected(
        "expected", Foam::dimensionSet(0, 0, -1, 0, 0),
                    3.199680897991421460e-21);

    Foam::coalescenceKernels::CoulaloglouTavlaridesCImpl kernel(
                C1, C2, alpha, sigma);

    auto result = kernel.S(xi1, xi2, epsilon, rho_d, nud);

    ASSERT_EQ(expected.value(), result.value());
    ASSERT_EQ(expected.dimensions(), result.dimensions());
}
