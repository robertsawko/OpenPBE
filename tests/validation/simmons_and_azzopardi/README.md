# Simmons and Azzopardi case

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
