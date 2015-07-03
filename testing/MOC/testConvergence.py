from numpy import (
    arange, array, exp, genfromtxt, linspace, nonzero,
    sqrt, piecewise, trapz, zeros)
from numpy.testing import assert_almost_equal, assert_array_less
from os import listdir
import re
from scipy.special import gamma


def scott_total_number_solution3(t, C=1, N0=1):
    T = C * N0 * t
    return 2.0 * N0 / (T + 2.0)


def scott_pbe_solution3(xi, t, C=1, N0=1, xi0=1):
    T = C * N0 * t
    x = xi / xi0
    phi3 = sum([
        (x * 2)**(2 * (k + 1)) / gamma(2 * (k + 1)) * (T / (T + 2))**k
        for k in range(100)
    ]) * 4.0 * exp(- 2 * x) / (x * (T + 2)**2)

    return N0 / xi0 * phi3


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


def test_breakup_convergence():
    from setupCases import breakup_cases

    L2_pbe_errors = zeros(len(breakup_cases))
    T = 10

    for k, bc in enumerate(breakup_cases):
        # Loading only number functions for each class
        path = "{0}{1}/postProcessing/probes/0/".format(bc.name, bc.nr_classes)
        classes = listdir(path)

        data = dict(
            (
                int(re.findall('\d+', c)[0]),
                genfromtxt(path + c)[:, 1]
            ) for c in classes
        )

        xi_n = linspace(bc.dv, bc.nr_classes * bc.dv, bc.nr_classes)
        Nsimulation = array([data[i][-1] for i in range(bc.nr_classes)])
        L2_pbe_errors[k] = L2_relative_error(
            xi_n,
            Nsimulation / bc.dv,
            zm_pure_breakup_pbe_solution(xi_n, T, bc.vmax)
        )

    check_error_convergence(L2_pbe_errors)


def test_coalescence_convergence():
    from setupCases import coalescence_cases

    L2_pbe_errors = zeros(len(coalescence_cases))
    T = 1

    for k, cc in enumerate(coalescence_cases):
        # Loading only number functions for each class
        path = "{0}{1}/postProcessing/probes/0/".format(cc.name, cc.nr_classes)
        classes = listdir(path)

        data = dict(
            (
                int(re.findall('\d+', c)[0]),
                genfromtxt(path + c)[:, 1]
            ) for c in classes
        )

        xi_n = linspace(cc.dv, cc.nr_classes * cc.dv, cc.nr_classes)
        Nsimulation = array([data[i][-1] for i in range(cc.nr_classes)])
        L2_pbe_errors[k] = L2_relative_error(
            xi_n,
            Nsimulation / cc.dv,
            scott_pbe_solution3(xi_n, T, cc.C, 2.0 * cc.v0, cc.N0)
        )

    check_error_convergence(L2_pbe_errors)
