#include <gtest/gtest.h>
#include <cmath>
#include "Integrator.H"


TEST(Integration, exponential_inverse){

    auto f = [](double x){return std::exp(-x);};

    auto r = integrate(f, 0.);

    ASSERT_DOUBLE_EQ(r, 1.);
}

