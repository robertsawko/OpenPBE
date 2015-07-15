#include <gtest/gtest.h>
// #include <cmath>
#include "scalar.H"
#include "../src/PBESystem/QMOM.H"
#include <Eigen/Dense>

/*
Unit tests for Wheeler moment inversion algorithm.
Based on examples 3.1 and 3.2 from:
Daniele L. Marchisio, Rodney O. Fox "Computational Models for Polydisperse
Particulate and Multiphase Systems", Cambridge University Press 2013
*/

//using Foam::sqrt;
//using Foam::pow;
//using Foam::scalar;
using Eigen::VectorXd;


TEST(qmom, wheeler_test){

    VectorXd moments(8);
    moments << 1, 5, 26, 140, 778, 4450, 26140, 157400;
    // Analytical solution
    double w1 = 1 / (2 * sqrt(6) + pow(2, 1.5) * sqrt(3) + 12);
    double w2 = 1 / (12 - 2 * sqrt(6) - pow(2, 1.5) * sqrt(3));
    VectorXd weights(4), abcissas(4);
    weights <<  w1, w2, w2, w1;
    abcissas << 2.6656, 4.2580, 5.7420, 7.3344;

    auto p = wheeler_inversion(moments);
    
    ASSERT_DOUBLE_EQ((p.abcissas - abcissas).norm(), 0.);
    ASSERT_DOUBLE_EQ((p.weights - weights).norm(), 0.);

}
