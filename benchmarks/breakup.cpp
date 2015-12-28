#include <benchmark/benchmark.h>
#include "breakupKernels/binaryBreakup.H"
#include "breakupKernels/noBreakup.H"
#include "breakupKernels/CoulaloglouTavlarides.H"

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

static void CoulaloglouTavlaridesBreakup(benchmark::State& state) {
    Foam::scalar C1(0.0152), C2(0.0678),
                 sigma(0.04282), alpha(0.050000000000000003);

    Foam::scalar epsilon(0.12935361098784559), rho_d(972.0),
           xi(9.738310380055584507e-11);

    Foam::breakupKernels::CoulaloglouTavlaridesImp kernel(C1, C2, alpha, sigma);

    while (state.KeepRunning())
        kernel.S(xi, rho_d, epsilon);
}

BENCHMARK(CoulaloglouTavlaridesBreakup);
