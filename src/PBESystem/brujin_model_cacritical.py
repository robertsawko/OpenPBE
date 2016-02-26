'''
The function for calculating critical capillary numbers.
Implementation based on description published in:

    Hill, D.P. "The Computer Simulation fo Dispersed Two-phase
    flows" PhD thesis, Imperial College, p 121

For viscous turbulent flows the flow parameter eta is 0.2
'''
from numpy import exp, log, sqrt

table32 = {
    0.0: (0.9353264, 0.5570426, -0.0000215, -0.0019007, -0.7542608),
    0.2: (0.3273221, 0.4726207, -0.0091768, -0.2481238, -1.3686552),
    0.4: (0.2858330, 0.3795495, -0.0171332, -0.2919881, -1.9388083),
    0.6: (0.2443439, 0.2864782, -0.0250896, -0.3358528, -2.5089614),
    0.8: (0.3102020, 0.2732713, -0.3524792, -0.4351167, -3.5251310),
    1.0: (0.4653474, 0.4538783, -0.5259108, -0.6653536, -5.2591139)}


def ca_critical(viscosity_ratio, eta=0.2):
    '''
    @param: viscosity_ratio molecular viscosity ratio
    '''
    c1, c2, c3, c4, c5 = table32[eta]
    logCa = sqrt(
        (c2 / c3 * log(viscosity_ratio) + c5 / c3)**2 -
        c1 / c3 * log(viscosity_ratio)**2 -
        2 * c4 / c3 * log(viscosity_ratio) -
        1 / c3) - c2 / c3 * log(viscosity_ratio) - c5 / c3 - 1
    return exp(logCa)
