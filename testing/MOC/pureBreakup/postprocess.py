import numpy as np
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt

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


def total_number_graph(N=10, resolution=300, deltaX=0.01):
    data = dict(
        (
            n,
            np.genfromtxt("postProcessing/probes/0/n{0:d}".format(n))[:, 1]
        ) for n in range(N))

    time = np.genfromtxt("postProcessing/probes/0/n0")[:, 0]
    N = sum(data.values())

    #set_plt_params()
    fig = plt.figure()
    ax = fig.gca()

    plt.xlabel("Time [s]")
    plt.ylabel("$N(t)/N(0)$")
    ax.loglog(time, N/N[0], label="Simulation")
    ax.legend(loc='lower right', shadow=True)
    fig.savefig(
        "total_number.pdf",
        bbox_inches='tight', dpi=300)
    plt.close()

total_number_graph()

#time_averaged_jointpdf(
    #N0=4000, N=10, resolution=0.02, a=3,
    #x1=np.array([0, 0, 0]),
    #x2=np.array([6, 0, 0]))
#time_averaged_pdf_product(
    #N0=4000, N=10, resolution=0.01, a=3,
    #x1=np.array([0, 0, 0]),
    #x2=np.array([6, 0, 0]))
