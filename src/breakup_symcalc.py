"""
Express the log-normal integral expressions as normal expressions. The latter
can be substituted then to CDF expressions.

Substitution is of the form:
x = log(d)
dx = 1/d dd

"""
from sympy import symbols, Rational, simplify, expand, init_printing
from sympy.assumptions import assuming, Q
from sympy.abc import mu, sigma, gamma

init_printing()

x = symbols('x', real=True)
with assuming(Q.integer(gamma), Q.positive(sigma), Q.real(mu)):
    expression = 2 * sigma**2 * x * (gamma - Rational(3, 2)) - (x - mu)**2

    # We are looking for a an expression of type (x - a)**2 -/+ c. We know that
    # a is
    a = mu + (gamma - Rational(3, 2)) * sigma**2
    a2 = a**2
    print(simplify(expand(a2) - mu**2))
