#ifndef ASTRO_HARSHDAS_GROWTH_FACTOR_H
#define ASTRO_HARSHDAS_GROWTH_FACTOR_H

#include "astro/cosmology/cosmologyBase.h"
#include "astro/cosmology/abstractGrowthFactor.h"
#include "astro/cosmology/gEffective.h"

namespace astro
{
  /**
   * \ingroup cosmologicalModel
   *
   * \class harshdasGrowthFactor
   *
   * \brief Class defining the standard cosmological growth factor
   *
   * This class defines the standard cosmological growth factor, i.e. for
   * laws of gravity which scale with wave number as \f$k^{-2}\f$. The class
   * is derived from the abstract class growthFactor.
   *
   * \author Matthias Bartelmann, Heidelberg University
   */
  class harshdasGrowthFactor: public abstractGrowthFactor
  {
  protected:
    cosmologyBase * cosmological_model;
    gEffective * g_effective;
    bool g_eff_internal;
    double k_internal;

    /**
     * Provides the differential equation to be solved
     */
    int differential_equation (double x, const double * y, double * dy);
  public:
    /**
     * Constructor initializing the growth factor with a pointer to a
     * cosmological model and optionally an effective gravitational
     * constant.
     * \param cosmological_model_in pointer to input cosmological model
     * \param g_effective_in pointer to input model for effective
     * gravitational constant
     */
    harshdasGrowthFactor
      (cosmologyBase * cosmological_model_in,
       gEffective * g_effective_in = NULL);

    /**
     * Destructor
     */
    ~harshdasGrowthFactor ();

    /**
     * Function returning the growth factor as a function of the scale
     * factor \c a, divided by the scale factor,
     * \f[
     *   G_+(a,k) = \frac{D_+(a,k)}{a}\;.
     * \f]
     * \param a scale factor
     * \param k wave number
     */
    double operator () (double a, double k = 1.0);

    /**
     * Returns the effective gravitational constant at scale factor \c a
     * \param a scale factor
     */
    double gravitationalConstant (const double a, double k);

    /**
     * Dummy function to comply with abstract base class. Always gives error
     */
    double gravitationalConstant (const double a);
  };
}

#endif
