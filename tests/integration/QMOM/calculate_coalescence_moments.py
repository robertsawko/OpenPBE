from sympy import symbols, sympify, integrate, summation, exp, gamma, oo
from numpy import linspace, zeros
from fractions import Fraction
import pickle

C = 0.1
tmax = 1
N = 100
M = 5
t = linspace(0, tmax, N)
m = zeros((N, M))
N0, v0 = symbols('N0 v_0', real=True, positive=True)
v, k = symbols('v k', real=True, positive=True)
x = v / v0

for i in range(N):
    T = Fraction(C * t[i] * 2).limit_denominator()
    for j in range(M):
        integral = integrate(
            v**j *
            N0 / v0 *
            4 * exp(-2 * x) *
            (T + 2)**(-2) *
            x**(-1) *
            (2 * x)**(2 * (k + 1)) / gamma(2*k+2),
            (v, 0, oo))
        m[i, j] = summation(
                integral * sympify(T / (T + 2))**k,
                (k, 0, oo)).subs(v0, 0.5).subs(N0, 2)


with open('m_analytical.pickle', 'wb') as f:
    pickle.dump(m, f)
