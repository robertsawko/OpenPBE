#include <benchmark/benchmark.h>
#include "../src/PBESystem/QMOM.H"
#include <Eigen/Dense>

static void BM_Wheeler(benchmark::State& state) {
    VectorXd moments(8);
    moments << 1, 5, 26, 140, 778, 4450, 26140, 157400;
    while (state.KeepRunning())
        wheeler_inversion(moments);
}

BENCHMARK(BM_Wheeler);

BENCHMARK_MAIN();
