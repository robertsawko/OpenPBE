#include <gtest/gtest.h>
#include <cmath>


TEST(Integration, exponential_inverse){

    auto f = [](double x){return std::exp(-x);};

    auto r = integrate(f);

    ASSERT_DOUBLE_EQ(r, 1.);
}

