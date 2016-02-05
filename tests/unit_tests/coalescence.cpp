#include <gtest/gtest.h>
#include "coalescenceKernels/constantCoalescence.H"
#include "coalescenceKernels/noCoalescence.H"
#include "coalescenceKernels/CoulaloglouTavlarides.H"

TEST(coalescence, constant) {
    Foam::dimensionedScalar xi1("xi1", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::dimensionedScalar xi2("xi2", Foam::dimensionSet(0, 3, 0, 0, 0), 8);

    Foam::scalar c = 5;

    Foam::coalescenceKernels::constantCoalescenceImpl kernel(c);

    Foam::dimensionedScalar expected(
        "expected", Foam::dimensionSet(0, 0, -1, 0, 0), c);

    auto result = kernel.S();

    ASSERT_EQ(expected.value(), result);
    // ASSERT_EQ(expected.dimensions(), result.dimensions());
}

TEST(coalescence, none) {
    Foam::dimensionedScalar xi1("xi1", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::dimensionedScalar xi2("xi2", Foam::dimensionSet(0, 3, 0, 0, 0), 8);

    Foam::coalescenceKernels::noCoalescenceImpl kernel;

    Foam::dimensionedScalar expected(
        "expected", Foam::dimensionSet(0, 0, -1, 0, 0), 0);

    auto result = kernel.S();

    ASSERT_EQ(expected.value(), result);
    // ASSERT_EQ(expected.dimensions(), result.dimensions());
}

TEST(coalescence, CoulaloglouTavlarides) {
    Foam::scalar C1(1.06e-12), C2(51300000000000.0), sigma(0.04282),
        alpha(0.050000000000000003);

    Foam::dimensionedScalar epsilon(
        "epsilon", Foam::dimensionSet(0, 2, -3, 0, 0), 0.12935361098784559);
    Foam::dimensionedScalar rho_c(
        "rho_c", Foam::dimensionSet(1, -3, 0, 0, 0), 1000);

    Foam::dimensionedScalar nuc(
        "nu_c", Foam::dimensionSet(0, 2, -1, 0, 0), 1e-6);

    Foam::dimensionedScalar xi1(
        "xi1", Foam::dimensionSet(0, 3, 0, 0, 0), 5.071423503403169858e-10);

    Foam::dimensionedScalar xi2(
        "xi2", Foam::dimensionSet(0, 3, 0, 0, 0), 2.885425297794245861e-11);

    Foam::dimensionedScalar Vcell(
        "xi2", Foam::dimensionSet(0, 3, 0, 0, 0), 0.012);

    Foam::dimensionedScalar expected("expected",
                                     Foam::dimensionSet(0, 0, -1, 0, 0),
                                     3.199680897991421460e-21);

    Foam::coalescenceKernels::CoulaloglouTavlaridesCImpl kernel(
        C1, C2, alpha, sigma);

    auto result = kernel.S(xi1.value(),
                           xi2.value(),
                           epsilon.value(),
                           rho_c.value(),
                           nuc.value(),
                           Vcell.value());

    ASSERT_FLOAT_EQ(expected.value(), result);
    // for (int i = 0; i < 7; ++i)
    //    ASSERT_NEAR(expected.dimensions()[i], result.dimensions()[i], 1e-9);
}
