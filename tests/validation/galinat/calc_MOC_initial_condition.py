from numpy import genfromtxt, pi, linspace, cumsum, clip
from scipy.interpolate import interp1d
from scipy.integrate import trapz


def calculate_MOC_initial_condition(number_of_classes, U=0.177, D=0.03):
    data = genfromtxt("paper_data/re6400/upstream.csv")

    # Data provided are diameters in mm. We need to covert to volume
    mm2m = 0.001
    v_data = pi / 6 * (data[:, 0] * mm2m)**3
    p_data = data[:, 1]

    p_interpolated = interp1d(
        v_data, p_data, kind='nearest', fill_value=0, bounds_error=False)

    # From Galinat we know that phase volumetric fraction is 1.7%-3%
    alpha = (0.03 + 0.017) / 2
    A = pi / 4 * D**2
    Vc = 1.0 * A  # Assume unit height column

    # Estimate the total number from total concentration
    v = linspace(0, 1.2 * max(v_data), number_of_classes + 1)
    dv = v[1]
    p_empirical = p_interpolated(v) / trapz(p_interpolated(v), x=v)
    cdf_empirical = cumsum(p_empirical) * dv  # Cumulative distribution

    meanVd = trapz(v * p_empirical, x=v)  # Dispersed phase mean volume
    N0 = alpha * Vc / meanVd  # Number of drops in 1m column

    # print("Total number of drops in unit column {0:0.0f}".format(N0))
    # print("Particle influx will be: {0:0.2f} drops/s".format(U * N0))
    # print(U * N0 * meanVd / (U * A))

    return dv, clip(
        N0 * (cdf_empirical[1:] - cdf_empirical[:-1]),
        a_min=1e-4, a_max=None
    )
