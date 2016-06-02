"""
Provide a volume distribution consistent with the log-normal diameter
distribution.

Sztencel and Jakubowski "Wstep do Rachunku Prawdopodobienstwa" section 5.5.
2000 Warszawa

If X has a continuous PDF f on (a,b) and function phi:(a, b) -> R is of class
C1 then Y=phi(X) has a PDF given by:

    g(y) = f(h(y)) |h'(y)| 1_phi(I)(y)

where, h(s) = phi^{-1}(s) and I=(a, b)
"""
from sympy import symbols, diff, pi, sqrt, exp, log, Rational


def construct_pdfs(mu, sigma):
    v = symbols('v', positive=True)
    d = symbols('d', positive=True)

    def diameter_pdf(d):
        return 1 / (d * sigma * sqrt(2 * pi)) * exp(
            - (log(d) - mu)**2 / (2 * sigma**2))

    h = (6 / pi * v)**Rational(1, 3)
    hprime = diff(h)
    return v, diameter_pdf(h) * abs(hprime), d, diameter_pdf(d)
