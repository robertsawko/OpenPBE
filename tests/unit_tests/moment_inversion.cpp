#include <gtest/gtest.h>
#include "scalar.H"
#include "MomentInversion.H"
#include <Eigen/Dense>

/*
Unit tests for Wheeler moment inversion algorithm.
Based on examples 3.1 and 3.2 from:
Daniele L. Marchisio, Rodney O. Fox "Computational Models for Polydisperse
Particulate and Multiphase Systems", Cambridge University Press 2013
*/

using Eigen::VectorXd;

TEST(wheeler_inversion, example31){

    VectorXd moments(8);
    moments << 1, 5, 26, 140, 778, 4450, 26140, 157400;
    // Analytical solution
    double w1 = 1 / (4 * sqrt(6) + 12);
    double w2 = 1 / (12 - 4 * sqrt(6));
    VectorXd weights(4), abcissas(4);
    weights <<  w1, w2, w2, w1;
    abcissas << 2.6656, 4.2580, 5.7420, 7.3344;

    auto result = wheeler_inversion(moments);
    
    ASSERT_NEAR((result.abcissas - abcissas).norm(), 0, 1e-4);
    ASSERT_NEAR((result.weights - weights).norm(), 0, 1e-4);
}

TEST(wheeler_inversion, example32){

    VectorXd moments(8);
    moments << 1, 0, 1, 0, 3, 0, 15, 0;
    // Analytical solution
    double xi1 = sqrt(sqrt(6) + 3);
    double xi2 = sqrt(3 - sqrt(6));
    double w1 = 1 / (4 * sqrt(6) + 12); 
    double w2 = 1 / (12 - 4 * sqrt(6)); 
    VectorXd weights(4), abcissas(4);
    weights <<  w1, w2, w2, w1;
    abcissas << -xi1, -xi2, xi2, xi1;

    auto result = wheeler_inversion(moments);
    
    ASSERT_NEAR((result.abcissas - abcissas).norm(), 0, 1e-4);
    ASSERT_NEAR((result.weights - weights).norm(), 0, 1e-4);
}
