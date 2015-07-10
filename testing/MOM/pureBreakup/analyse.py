from numpy import genfromtxt, array, dot, column_stack, vstack, append, linspace
from setupIC import class_number, v, dv, m0, prob
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt

path = "./MOC{0}/postProcessing/probes/0/".format(class_number)

MOCdata = genfromtxt(path + 'n0')

for k in range(1, class_number):
    MOCdata = column_stack(
        (MOCdata, genfromtxt(path + 'n{0}'.format(k))[:, 1])
    )

path = "./postProcessing/probes/0/"

MOMdata = genfromtxt(path + 'm0', usecols=1)
for k in range(1, 3):
    MOMdata = column_stack(
        (MOMdata, genfromtxt(path + 'm{0}'.format(k))[:, 1])
    )

NDFinit = m0 * (prob.cdf(v) - prob.cdf(v - dv)) / dv
NDFend = MOCdata[-1][1:] / dv
time = append([0.0], MOCdata[:, 0])
m = MOMdata[-1, :]
k = m[1]**2 / (m[2] * m[0] - m[1]**2)
theta = m[0] * m[1] / (m[0] * m[2] - m[1]**2)

from scipy.stats import gamma
prob_MOM = gamma(a=k, scale=theta)

NDFMOM = m[0] * (prob_MOM.pdf(v) - prob_MOM.cdf(v - dv)) / dv

fig = plt.figure()
ax = fig.gca()
plt.xlabel("Drop volume")
plt.ylabel("Number density")
ax.plot(v, NDFinit)
ax.plot(v, NDFend)
ax.plot(v, NDFMOM)
fig.savefig("pbe.pdf", bbox_inches='tight')


def calculate_moment(N, k=0):
    return dot(v**k, N) * dv

MOC_moments = array([calculate_moment(NDFinit, k) for k in range(3)])
for n in range(MOCdata.shape[0]):
    ndf = MOCdata[n][1:] / dv
    MOC_moments = vstack(
        (
            MOC_moments,
            array([calculate_moment(ndf, k) for k in range(3)])
        )
    )



for k in range(3):
    fig = plt.figure()
    ax = fig.gca()
    plt.xlabel('Time [s]')
    # plt.ylabel('Moment {0}'.format(k))
    ax.plot(time, MOC_moments[:, k])
    ax.plot(time[1:], MOMdata[:, k])
    fig.savefig('moment{0}.pdf'.format(k), bbox_inches='tight')


