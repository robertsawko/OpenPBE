#include <benchmark/benchmark.h>
#include "coalescenceKernels/constantCoalescence.H"
#include "coalescenceKernels/noCoalescence.H"

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
