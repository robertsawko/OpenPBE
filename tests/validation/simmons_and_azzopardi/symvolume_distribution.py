"""
Provide a volume distribution consistent with the log-normal diameter
distribution.

Sztencel and Jakubowski "Wstęp do Rachunku Prawdopodobieństwa" 5.5.

If X has a continuous PDF f on (a,b) and function phi:(a, b) -> R is of class
C1 then Y=phi(X) has a PDF given by:

    g(y) = f(h(y)) |h'(y)| 1_phi(I)(y)

where, h(s) = phi^{-1}(s) and I=(a, b)
"""
from sympy import symbols, diff, pi, sqrt, exp, log, Rational, Function
from sympy.utilities.lambdify import implemented_function
from sympy import pprint

v = symbols('v', positive=True)
sigma = 0.23482069474132783
mu = -7.8028262260262604


def diameter_pdf(d):
    return 1 / (d * sigma * sqrt(2 * pi)) * exp(
        - (log(d) - mu)**2 / (2 * sigma**2))

h = (6 / pi * v)**Rational(1, 3)
hprime = diff(h)
volume_pdf = diameter_pdf(h) * abs(hprime)
pprint(volume_pdf)
