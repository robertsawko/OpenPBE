from sympy import init_session
from fractions import Fraction

init_session()
N0, v0 = symbols("N_0 v_0", real=True, positive=True)
v = symbols('v', real=True, positive=True)
k = symbols('k', integer=True, positive=True)
x = v / v0

expressions = dict()
T = sympify(Fraction(2, 10))

for m in range(5):
    integral = integrate(
        v**m *
        N0 / v0 *
        4 * exp(-2 * x) *
        (T + 2)**(-2) *
        x**(-1) *
        (2*x)**(2*(k+1)) / gamma(2*k+2),
        (v, 0, oo))

    expressions[m] = summation(integral * (T / (T + 2))**(k), (k, 0, oo))

display(expressions)
