#ifndef INTEGRATOR_H
#define INTEGRATOR_H
#include <gsl/gsl_integration.h>

template <typename F>
double toIntegrate(double x, void* params)
{
    auto f = static_cast<F*>(params);
    return (*f)(x);
}

template <typename F>
double integrate(F f, double lower_bound)
{
    auto integrationWorkspace = gsl_integration_workspace_alloc(1000);
    gsl_function gslFunction;
    gslFunction.function = &toIntegrate<F>;
    gslFunction.params = &f;

    double result, error;
    gsl_integration_qagiu(&gslFunction, lower_bound, 1e-3, 0, 1000, integrationWorkspace,
                          &result, &error);

    gsl_integration_workspace_free(integrationWorkspace);

    return result;
}

#endif // INTEGRATOR_H

