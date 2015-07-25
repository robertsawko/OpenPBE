#include <benchmark/benchmark.h>
#include "../src/PBESystem/breakupKernels/binaryBreakup.H"

static void binaryBreakup(benchmark::State& state) {
    Foam::dimensionedScalar xi(
                "xi", Foam::dimensionSet(0, 3, 0, 0, 0), 6);

    Foam::breakupKernels::binaryBreakup kernel;

    while (state.KeepRunning())
        kernel.S(xi);
}

BENCHMARK(binaryBreakup);
