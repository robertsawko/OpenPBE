# Simmons and Azzopardi

Simmons, M. and Azzopardi, J. "Drop size distributions in
dispersed liquid--liquid pipe flow", International Journal of
Multiphase Flow, 2001

Pipe diameter is 0.063m

Horizontal case is selected

| Dir | Mixture superficial velocity [m/s] | Dispersed phase concentration [%]|
| --- |:----------------------------------:|:--------------------------------:|
| c01 | 2.551                              | 6.2                              |

## Inlet and initial conditions

The inlet is subivided into a ring area and the centre area. The ring area
occupies the area equal to the percentage of the volume fraction of the
dispersed phase.

### Velocity

Mixture velocity is specified for both inlet regions.

### Method of classes

Inlet and initial conditions:
 * log-normally distributed diameter,
 * mean diameter 500micron with 50micron std,
 * consistent recalculation for volume,
 * total number function used

At the moment we're using total number rather than number density. It seems
that that is related to the way we define source term and set up initalisation.
The current implementation of CT is for total number!

## Content
 * `symcalc.py` - gives the form of diameter and volume distributions
   consistent with log-normal diameter distribution assumption.
 * `compare_distributions.py` - a script comparing the distribution with the
   use of sampling.
 * `setupCases.py` - a script setting up the inlet and initial condition.
