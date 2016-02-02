from numpy import pi, exp, log, sqrt
from scipy.stats import lognorm
from sympy import lambdify
from numpy import histogram, linspace
from symvolume_distribution import volume_pdf, v
from matplotlib.pyplot import figure

# Below are the sample mean and standard deviation of the non-logorithmicied
# distribution
m = 4.2e-4
s = 1e-4
s2 = s**2
# Using the moment inversion found e.g.
# https://en.wikipedia.org/wiki/Log-normal_distribution#Notation
mu = log(m / sqrt(1 + s2 / (m**2)))
sigma = sqrt(log(1 + s2 / m**2))

def diameter_pdf(d):
    return 1 / (d * sigma * sqrt(2 * pi)) * exp(
        - (log(d) - mu)**2 / (2 * sigma**2))

f = lambdify(v, volume_pdf, "numpy")

N = 100000
d_samples = lognorm(s=sigma, scale=exp(mu)).rvs(N)
v_samples = pi / 6 * d_samples**3


def plot(samples, expression, max_val, name):
    numbers, edges = histogram(samples, bins=100, range=(0, max_val))
    width = edges[1] - edges[0]
    grid = linspace(edges[0], edges[-1], 1000)
    fig = figure()
    ax = fig.gca()
    ax.set_xlabel(name)
    ax.bar(
        edges[:-1], numbers / N / width, width=0.9 * width, alpha=0.3,
        label="Sampling")
    ax.plot(grid, expression(grid), 'g', linewidth=3, label="Expression")
    fig.savefig('{0:s}.png'.format(name))


plot(d_samples, diameter_pdf, max_val=1e-3, name="Diameter")
plot(v_samples, f, max_val=5e-10, name="Volume")
