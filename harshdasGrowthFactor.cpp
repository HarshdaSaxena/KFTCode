#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "harshdasGrowthFactor.h"
#include "gNewton.h"
#include "../utilities/deSolver.h"

astro::harshdasGrowthFactor::harshdasGrowthFactor
  (cosmologyBase * cosmological_model_in, astro::gEffective * g_effective_in):
cosmological_model (cosmological_model_in), g_effective (g_effective_in)
{
  if (g_effective_in == NULL)
  {
    g_eff_internal = true;
    g_effective = new astro::gNewton;
  }
  else
  {
    g_eff_internal = false;
    g_effective = g_effective_in;
  }

  k_internal = 1.0;
}

astro::harshdasGrowthFactor::~harshdasGrowthFactor ()
{
  if (g_eff_internal)
    delete g_effective;
}

int astro::harshdasGrowthFactor::differential_equation
  (double x, const double y[], double dy[])
{
  double c[2] =
   {3.0/x+cosmological_model->dexpansionFunction_da (x)/
    cosmological_model->expansionFunction (x),
    1.5*g_effective->operator () (x,k_internal)*cosmological_model->omegaMatter (x)/x/x};
  dy[0] = y[1];
  dy[1] = -y[1]*c[0]+y[0]*c[1];
  return GSL_SUCCESS;
}

double astro::harshdasGrowthFactor::operator () (double a, double k)
{
  double a_equality = cosmological_model->getEqualityScale ();
  double n_growth = cosmological_model->getGrowthExponent ();
  if (a < a_equality)
    return pow (a_equality, n_growth-1.0);

  k_internal = k;

  std::function<int (double, const double[], double[])> diffeq = std::bind
    (&astro::harshdasGrowthFactor::differential_equation, this,
     std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);
  astro::deSolver solve (diffeq, 2);
  double y[2] = {a_equality, n_growth*pow (a_equality, n_growth-1.0)};
  std::vector<double> result = solve (y, a_equality, a);
  return result[0]/a;
}

double astro::harshdasGrowthFactor::gravitationalConstant (double a, double k)
{ return g_effective->operator () (a, k); }
