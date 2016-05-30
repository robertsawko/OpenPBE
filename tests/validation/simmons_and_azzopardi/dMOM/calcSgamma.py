from numpy import log, sqrt, pi, mean
from numpy.random import lognormal

dmean = 5.0e-4  # 500 microns
dstd = 5.0e-5  # std of 10% from the mean
dvar = dstd**2
mu = log(dmean / sqrt(1 + dvar / (dmean**2)))
sigma = sqrt(log(1 + dvar / dmean**2))

# Monte carlo calculation for mean volume
ds = lognormal(mean=mu, sigma=sigma, size=100000)
vmean = mean(pi / 6 * ds**3)
# Calculate total number of drops
dpipe = 0.063  # [m]
alpha = 0.062  # [1]
S0 = alpha / vmean
S1 = S0 * dmean
S2 = S0 * (dvar + dmean**2)
print(S0, S1, S2)
print(mu, sigma)
print(S2 / S0)
