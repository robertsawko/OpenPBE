#include <gtest/gtest.h>
#include "coalescenceKernels/constantCoalescence.H"
#include "coalescenceKernels/noCoalescence.H"

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
