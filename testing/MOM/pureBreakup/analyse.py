from numpy import genfromtxt, array, dot
from os import listdir
import re
from setupIC import class_number, v, dv, m0, prob
import matplotlib as mpl
mpl.use("agg")
import matplotlib.pyplot as plt

path = "./MOC{0}/postProcessing/probes/0/".format(class_number)
classes = listdir(path)

data = dict(
    (
        int(re.findall('\d+', c)[0]),
        genfromtxt(path + c)[:, 1]
    ) for c in classes
)

Ninit = m0 * (prob.cdf(v) - prob.cdf(v - dv)) / dv
Nend = array([data[i][-1] for i in range(class_number)]) / dv


def calculate_moment(N, k=0):
    return dot(v**k, N) * dv

fig = plt.figure()
ax = fig.gca()
plt.xlabel("Drop volume")
plt.ylabel("Number density")
ax.plot(v, Ninit)
ax.plot(v, Nend)
fig.savefig("pbe.pdf", bbox_inches='tight')

for k in range(3):
    print(k)
    print(calculate_moment(Ninit, k))
    print(calculate_moment(Nend, k))
