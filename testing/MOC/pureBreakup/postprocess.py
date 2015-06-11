import numpy as np
import matplotlib as mpl
mpl.use("agg")

import matplotlib.pyplot as plt
from itertools import cycle
from os import listdir
import re

plt.style.use('ggplot')
plt.ioff()


def set_plt_params(
        relative_fig_width=1.0, landscape=True, page_width=307.3, rescale_h=1):

    fig_width_pt = page_width * relative_fig_width
    inches_per_pt = 1.0 / 72.27               # Convert pt to inch
    golden_mean = (np.sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches

    if landscape:
        fig_height = fig_width * golden_mean  # height in inches
    else:
        fig_height = fig_width / golden_mean  # height in inches

    fig_height = fig_height * rescale_h
    fig_size = [fig_width, fig_height]
    params = {
        'font.family': 'serif',
        'axes.labelsize': 7,
        'xtick.labelsize': 5,
        'ytick.labelsize': 5,
        'axes.labelcolor': 'black',
        'ytick.color': 'black',
        'xtick.color': 'black',
        'legend.handlelength': 4,
        'legend.fontsize': 7,
        # 'lines.markersize': 3,
        # 'xtick.labelsize': 7,
        # 'ytick.labelsize': 7,
        'text.usetex': True,
        'text.latex.unicode': True,
        'figure.figsize': fig_size,
        'pgf.texsystem': "xelatex",
        'pgf.rcfonts': False,
    }

    plt.rcParams.update(params)


def zm_pure_breakup_total_number_solution(x, t, l):
    """
    This is simply an integral of Ziff and McGrady
    """
    return np.exp(-t * l**2) \
        + np.trapz(2.0 * t * l * np.exp(-t * x**2), x=x)


def zm_pure_breakup_pbe_solution(x, t, l):
    """
    This is based on Equation 25 from Ziff and McGrady
    """
    return np.piecewise(
        x,
        [x < l, x == l, x > l],
        [
            lambda x: 2.0 * t * l * np.exp(-t * x**2),
            lambda x: np.exp(-t * x**2),
            lambda x: 0.0
        ]

    )


from parameterizedVariation import nr_classes, l


total_fig = plt.figure()
total_ax = total_fig.gca()
plt.xlabel("Time [s]")
plt.ylabel("$N(t)/N(0)$")
pbe_fig = plt.figure()
pbe_ax = pbe_fig.gca()
plt.xlabel("Drop volume")
plt.ylabel("Number density")
linestyles = cycle(['-', '--', ':'])
ts = [5, 10]

# Will use opacity but with consistent colors for each time step
cc = pbe_ax._get_lines.color_cycle
colors = dict(
    [(t, next(cc)) for t in ts]
)

for n in nr_classes:
    # Loading only number functions for each class
    path = "testCase{0}/postProcessing/probes/0/".format(n)
    classes = listdir(path)

    data = dict(
        (
            int(re.findall('\d+', c)[0]),
            np.genfromtxt(path + c)[:, 1]
        ) for c in classes
    )

    time = np.genfromtxt(path + "n0")[:, 0]

    Ntotal = sum(data.values())
    total_ax.loglog(time, Ntotal / Ntotal[0], label="MOC N={0}".format(n))

    deltaX = 1.0 / n
    xi_n = np.linspace(deltaX, n * deltaX, n)

    for t in ts:
        ind = np.nonzero(time == t)[0]
        Nsimulation = np.array([data[i][ind] for i in range(n)])
        pbe_ax.semilogy(
            xi_n, Nsimulation / (l / n), "-",
            color=colors[t],
            alpha=n / 40.0,
            label="MOC N={0} for $t={1}$".format(n, t)
        )

xi_a = np.linspace(0, l, 100, endpoint=False)
N_analytical = dict(
    (t, zm_pure_breakup_pbe_solution(xi_a, t, l))
    for t in ts
)
for t in ts:
    pbe_ax.semilogy(
        xi_a, N_analytical[t], "k",
        label="Ziff and McGrady $t={0}$".format(t),
        linestyle=next(linestyles)
    )

pbe_ax.legend(loc='lower left', shadow=True)

# Evalute the analytical solution for total number
xi = np.linspace(0, l, 100)
N_analytical = np.zeros(time.shape)
for i, t in enumerate(time):
    N_analytical[i] = zm_pure_breakup_total_number_solution(xi, t, l)
total_ax.loglog(
    time, N_analytical / N_analytical[0], "k", label="Ziff and McGrady")
total_ax.legend(loc='upper left', shadow=True)

total_fig.savefig("total_number.pdf", bbox_inches='tight')
pbe_fig.savefig("pbe.pdf", bbox_inches='tight')
plt.close()
