#include <benchmark/benchmark.h>
#include "breakupKernels/binaryBreakup.H"
#include "breakupKernels/noBreakup.H"

static void binaryBreakup(benchmark::State& state) {
    Foam::dimensionedScalar xi(
                "xi", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::breakupKernels::binaryBreakupImpl kernel;

    while (state.KeepRunning())
        kernel.S(xi);
}

BENCHMARK(binaryBreakup);

static void noBreakup(benchmark::State& state) {
    Foam::dimensionedScalar xi(
                "xi", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::breakupKernels::binaryBreakupImpl kernel;

    while (state.KeepRunning())
        kernel.S(xi);
}

BENCHMARK(noBreakup);

