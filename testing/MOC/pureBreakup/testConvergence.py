from numpy import (
    arange, array, exp, genfromtxt, linspace, nonzero,
    sqrt, piecewise, trapz, zeros)
import matplotlib as mpl
mpl.use("agg")
from numpy.testing import assert_almost_equal, assert_array_less
from os import listdir
import re


def zm_pure_breakup_total_number_solution(x, t, l):
    """
    This is simply an integral of Ziff and McGrady
    """
    return exp(-t * l**2) \
        + trapz(2.0 * t * l * exp(-t * x**2), x=x)


def zm_pure_breakup_pbe_solution(x, t, l):
    """
    This is based on Equation 25 from Ziff and McGrady
    """
    return piecewise(
        x,
        [x < l, x == l, x > l],
        [
            lambda x: 2.0 * t * l * exp(-t * x**2),
            lambda x: exp(-t * x**2),
            lambda x: 0.0
        ]

    )


def L2_relative_error(x, f, g):
    return sqrt(trapz((f - g)**2, x=x)) / sqrt(trapz(f**2, x=x))


def check_error_convergence(L2_errors):
    # Testing convergence
    for k in arange(1, len(L2_errors)):
        assert_array_less(L2_errors[k], L2_errors[k - 1])

    # Testing convergence this will equal to less than 1% error
    assert_almost_equal(L2_errors[-1], 0.0, decimal=1)


def test_convergence():
    from parameterizedVariation import nr_classes, l

    L2_pbe_errors = zeros(len(nr_classes))
    T = 10

    for k, N in enumerate(nr_classes):
        # Loading only number functions for each class
        path = "testCase{0}/postProcessing/probes/0/".format(N)
        classes = listdir(path)

        data = dict(
            (
                int(re.findall('\d+', c)[0]),
                genfromtxt(path + c)[:, 1]
            ) for c in classes
        )

        #time = genfromtxt(path + "n0")[:, 0]

        #Ntotal = sum(data.values())

        deltaX = 1.0 / N
        xi_n = linspace(deltaX, N * deltaX, N)
        Nsimulation = array([data[i][-1] for i in range(N)])
        L2_pbe_errors[k] = L2_relative_error(
            xi_n,
            Nsimulation / (l / N),
            zm_pure_breakup_pbe_solution(xi_n, T, l)
        )

    check_error_convergence(L2_pbe_errors)
