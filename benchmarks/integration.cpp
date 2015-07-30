#include <benchmark/benchmark.h>
#include "../src/PBESystem/Integrator.H"
#include <cmath>

static void gsl_integration(benchmark::State& state) {
    auto f = [](double x){return std::exp(-x);};

    while (state.KeepRunning())
        auto r = integrate(f, 0.);
}

BENCHMARK(gsl_integration);
