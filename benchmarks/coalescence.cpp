#include <benchmark/benchmark.h>
#include "coalescenceKernels/constantCoalescence.H"
#include "coalescenceKernels/noCoalescence.H"
#include "coalescenceKernels/CoulaloglouTavlarides.H"

static void constantCoalescence(benchmark::State& state) {
    Foam::scalar c = 5;

    Foam::coalescenceKernels::constantCoalescenceImpl kernel(c);

    while (state.KeepRunning())
        kernel.S();
}

BENCHMARK(constantCoalescence);

static void noCoalescence(benchmark::State& state) {

    Foam::coalescenceKernels::noCoalescenceImpl kernel;

    while (state.KeepRunning())
        kernel.S();
}

BENCHMARK(noCoalescence);

static void CoulaloglouTavlaridesCoalescence(benchmark::State& state) {
    Foam::scalar C1(0.0152), C2(0.0678),
                 sigma(0.04282), alpha(0.050000000000000003);

    Foam::scalar epsilon(0.12935361098784559), rho_d(972.0), nud(0.001),
        xi1(5.071423503403169858e-10), xi2(2.885425297794245861e-11),
        Vcell(0.012);

    Foam::coalescenceKernels::CoulaloglouTavlaridesCImpl kernel(
                C1, C2, alpha, sigma);

    while (state.KeepRunning())
        kernel.S(xi1, xi2, epsilon, rho_d, nud, Vcell);
}

BENCHMARK(CoulaloglouTavlaridesCoalescence);
