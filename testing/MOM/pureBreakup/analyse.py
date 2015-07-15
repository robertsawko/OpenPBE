from numpy import genfromtxt, array, dot, column_stack, vstack, append
from setupIC import class_number, v, dv, m0, prob
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt


def read_data_from_probes(path, variable_name, number_of_probes):
    data = genfromtxt(path + variable_name + '0')
    for k in range(1, number_of_probes):
        data = column_stack(
            (data, genfromtxt(
                "{0}{1}{2}".format(path, variable_name, k))[:, 1])
        )

    return data


MOCdata = read_data_from_probes(
    "./MOC{0}/postProcessing/probes/0/".format(class_number),
    'n', class_number)

MOMdata = read_data_from_probes(
    "./postProcessing/probes/0/", 'm', 3
)

NDFinit = m0 * (prob.cdf(v) - prob.cdf(v - dv)) / dv
NDFend = MOCdata[-1][1:] / dv
time = append([0.0], MOCdata[:, 0])


def calculate_moment(N, k=0):
    return dot(v**k, N) * dv


fig = plt.figure()
ax = fig.gca()
plt.xlabel("Drop volume")
plt.ylabel("Number density")
ax.plot(v, NDFinit)
ax.plot(v, NDFend)
fig.savefig("pbe.pdf", bbox_inches='tight')

MOC_moments = array([calculate_moment(NDFinit, k) for k in range(3)])
for n in range(MOCdata.shape[0]):
    ndf = MOCdata[n][1:] / dv
    MOC_moments = vstack(
        (
            MOC_moments,
            array([calculate_moment(ndf, k) for k in range(3)])
        )
    )

print(MOC_moments[0, :])
fig = plt.figure()
plt.xlabel('Time [s]')
plt.ylabel('Moment error [%]')
ax = fig.gca()
lines = []
for k in range(3):
    error = \
        abs(MOC_moments[1:, k] - MOMdata[:, k + 1]) \
        / MOMdata[:, k + 1] * 100
    l, = ax.plot(time[1:], error)
    lines.append(l)
ax.legend(lines, ["$m_0$", "$m_1$", "$m_2$"], loc='upper left')
fig.savefig('moment_errors.pdf'.format(k), bbox_inches='tight')
