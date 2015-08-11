from numpy import array, genfromtxt, linspace, nonzero, sqrt
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt
from itertools import cycle
from os import listdir
import re
import pickle

plt.style.use('ggplot')
plt.ioff()


def set_plt_params(
        relative_fig_width=1.0, landscape=True, page_width=307.3, rescale_h=1):

    fig_width_pt = page_width * relative_fig_width
    inches_per_pt = 1.0 / 72.27               # Convert pt to inch
    golden_mean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
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
        'lines.markersize': 4,
        # 'xtick.labelsize': 7,
        # 'ytick.labelsize': 7,
        'text.usetex': True,
        'text.latex.unicode': True,
        'figure.figsize': fig_size,
        'pgf.texsystem': "xelatex",
        'pgf.rcfonts': False,
    }

    plt.rcParams.update(params)


vmax = 1e1
N0 = 2
v0 = 0.5
C = 0.1
pkl_file = open('m_analytical.pickle', 'rb')

m_analytical = pickle.load(pkl_file)


linestyles = cycle(['-', '--', ':'])

#for m in range(5):

# Will use opacity but with consistent colors for each time step
#cc = pbe_ax._get_lines.color_cycle
#colors = dict(
    #[(t, next(cc)) for t in ts]
#)

set_plt_params()
data = dict()
for n in range(1, 5):
    # Loading only number functions for each class
    path = "pure_coal{0}/postProcessing/probes/0/".format(n)
    moment_files = listdir(path)

    moments = dict(
        (
            int(re.findall('\d+', mf)[0]),
            genfromtxt(path + mf)[:, 1]
        ) for mf in moment_files
    )
    data[n] = moments

time = genfromtxt(path + "m0")[:, 0]

ta = linspace(0, 1, m_analytical.shape[0])
lines = []
for m in range(4):
    fig = plt.figure()
    ax = fig.gca()
    plt.xlabel('Time [s]')
    plt.ylabel('Moment {0}'.format(m))
    l, = ax.plot(
        ta[1:], m_analytical[1:, m],
        linewidth=0, marker='o', markersize=2, markevery=3)
    lines.append(l)
    l, = ax.plot(time, data[m + 1][m], alpha=0.5)
    lines.append(l)
    fig.savefig("moment{0}.pdf".format(m), bbox_inches='tight')


figlegend = plt.figure(figsize=(0.5, 0.2))
figlegend.legend(lines[0:2], ["Analytical", "QMOM"])
figlegend.savefig("legend.pdf", bbox_inches='tight')
