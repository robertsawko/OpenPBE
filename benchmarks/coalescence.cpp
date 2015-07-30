#include <benchmark/benchmark.h>
#include "coalescenceKernels/constantCoalescence.H"
#include "coalescenceKernels/noCoalescence.H"
#include "coalescenceKernels/CoulaloglouTavlarides.H"

static void constantCoalescence(benchmark::State& state) {
    Foam::dimensionedScalar xi1(
        "xi1", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::dimensionedScalar xi2(
        "xi2", Foam::dimensionSet(0, 3, 0, 0, 0), 8);

    Foam::scalar c = 5;

    Foam::coalescenceKernels::constantCoalescenceImpl kernel(c);

    while (state.KeepRunning())
        kernel.S(xi1, xi2);
}

BENCHMARK(constantCoalescence);

static void noCoalescence(benchmark::State& state) {
    Foam::dimensionedScalar xi1(
        "xi1", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::dimensionedScalar xi2(
        "xi2", Foam::dimensionSet(0, 3, 0, 0, 0), 8);

    Foam::coalescenceKernels::noCoalescenceImpl kernel;

    while (state.KeepRunning())
        kernel.S(xi1, xi2);
}

BENCHMARK(noCoalescence);

static void CoulaloglouTavlaridesCoalescence(benchmark::State& state) {
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

    Foam::coalescenceKernels::CoulaloglouTavlaridesCImpl kernel(
                C1, C2, alpha, sigma);

    while (state.KeepRunning())
        kernel.S(xi1, xi2, epsilon, rho_d, nud);
}

BENCHMARK(CoulaloglouTavlaridesCoalescence);
