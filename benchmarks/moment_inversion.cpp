#include <benchmark/benchmark.h>
#include "../src/PBESystem/MomentInversion.H"
#include <Eigen/Dense>

static void Wheeler(benchmark::State& state) {
    Eigen::VectorXd moments(8);
    moments << 1, 5, 26, 140, 778, 4450, 26140, 157400;
    while (state.KeepRunning())
        wheeler_inversion(moments);
}

BENCHMARK(Wheeler);
