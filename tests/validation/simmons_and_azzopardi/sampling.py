"""
Sampling example that provides validation of symbolic calculation from symcalc.
In this example N diameters is selected with lognormal distribution and then
transformed into volumes. Subsequently plots are constructed with analytical
and empircal distributions.
"""
from numpy import pi, exp, log, sqrt
from scipy.stats import lognorm
from sympy import lambdify
from numpy import histogram, linspace
from symcalc import construct_pdfs
from matplotlib.pyplot import figure

# Below are the sample mean and standard deviation of the non-logorithmicied
# distribution
m = 5.0e-4
s = 5e-5  # This corresponds to std of 10% mean value
s2 = s**2  # Variance
# Using the moment inversion found e.g.
# https://en.wikipedia.org/wiki/Log-normal_distribution#Notation
mu = log(m / sqrt(1 + s2 / (m**2)))
sigma = sqrt(log(1 + s2 / m**2))

v, volume_pdf, d, diameter_pdf = construct_pdfs(mu, sigma)

f = lambdify(v, volume_pdf, "numpy")
g = lambdify(d, diameter_pdf, "numpy")

N = 100000
d_samples = lognorm(s=sigma, scale=exp(mu)).rvs(N)
v_samples = pi / 6 * d_samples**3


def plot(samples, expression, max_val, name):
    numbers, edges = histogram(samples, bins=100, range=(0, max_val))
    width = edges[1] - edges[0]
    grid = linspace(edges[1], edges[-1], 1000)
    fig = figure()
    ax = fig.gca()
    ax.set_xlabel(name)
    ax.bar(
        edges[:-1], numbers / N / width, width=0.9 * width, alpha=0.3,
        label="Sampling")
    ax.plot(grid, expression(grid), 'g', linewidth=3, label="Expression")
    ax.legend()
    fig.savefig('{0:s}.png'.format(name))


plot(d_samples, g, max_val=1e-3, name="Diameter")
plot(v_samples, f, max_val=3e-10, name="Volume")
