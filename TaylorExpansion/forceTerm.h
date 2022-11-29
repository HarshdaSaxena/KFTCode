#ifndef KFT_FORCE_TERM_H
#define KFT_FORCE_TERM_H

#include <astro/utilities/functionTable.h>
#include <astro/cosmology/cosmologyBase.h>
#include <astro/cosmology/abstractGrowthFactor.h>
#include <astro/cosmology/standardGrowthFactor.h>
#include <astro/cosmology/cosmicStructures.h>
#include <astro/cosmology/tophatFilter.h>
#include <astro/cosmology/gaussFilter.h>
#include <astro/cosmology/bardeenPowerSpectrum.h>
#include <astro/cosmology/smithEtalPowerSpectrum.h>
#include <astro/cosmology/gEffective.h>

#include <KFT/particleDynamics/particleDynamics.h>

class forceTerm
{
protected:
  double a_min, a_max, k_min, k_max;
  double D_plus_min, D_plus_max, t_max, D_plus_norm;
  double amplitude, sigma_1, viscosity, tau, k_scale, k_star;
  double a_save, k_save;
  double f_k_save, f_r_save;
  unsigned int n_table;
  int proptype;
  double a_norm;

  astro::cosmologyBase * cosmological_model;
  astro::abstractGrowthFactor * growth_model;
  astro::cosmicStructures * cosmic_structures;
  astro::tophatFilter * tophat_filter;
  astro::gaussFilter * gauss_filter;
  astro::bardeenPowerSpectrum * linear_spectrum;
  astro::smithEtalPowerSpectrum * nonlinear_spectrum;
  astro::gEffective * g_effective;

  KFT::particleDynamics * propagator;
  KFT::particleDynamics * propagator_h;

  astro::functionTable
    dDplus_da_tab,
    sigma_J_sq_table,
    sigma_J_prime_sq_table;
/*
    sigma_J_sq_dot_table,
    dsigma_J_sq_dA_table,
    v12_parallel_table,
    f12_parallel_table;
*/

  bool cos_model_internal, g_eff_internal;

  void init_tau ();
  void init_tables (double f_k = 1.0, double f_r = 1.0);
  void init_tables_zeldovich (double f_k = 1.0, double f_r = 1.0);

  double displacement (double t);
  double interaction_term_hamilton (double k, double a);
  double interaction_term_normal (double k, double a);
  double taylorfactor (double k, double a);

public:

  /**
   * Constructor initialisting the class.
   * \param growth_model_in pointer to input growth factor
   * \param cosmological_model_in pointer to input cosmological model
   * \param proptype_in integer deciding if hamilton or zeldovich propagators
   * should be used, 0 for hamilton, 1 for zeldovich
   * \param g_effective_in pointer to graviational constant
   * \param a_norm_in scale factor at which the the power spectrum is to be
   * normalized
   */
  forceTerm (astro::abstractGrowthFactor * growth_model_in = NULL,
             astro::cosmologyBase * cosmological_model_in = NULL,
             astro::gEffective * g_effective_in = NULL,
             double a_norm_in = 1.0, int proptype = 1);

  /**
   * Destructor
   */
  ~forceTerm();

  /**
   * Returns the kernel function to be integrated with the power spectrum to
   * obtain the mean force between particle pairs.
   * \f[
   *   J(y, y_0) =
   *     1+\frac{1-y^2-y_0^2}{4y}
   *     \ln\frac{y_0^2+(1+y)^2}{y_0^2+(1-y)^2}
   * \f]
   * \param y dimension-less wave number
   * \param y0 dimension-less Yukawa scale
   */
  double J (double y, double y0);
  double J_prime (double y, double y0, double y_star);

  /**
   * Returns the moment
   * \f[
   *   \sigma_J^2 =
   *     \frac{k^3}{(2\pi)^2}\int_0^\infty dy\,y^2\,P(ky)\,J(y,y_0)
   * \f]
   * of the initial power spectrum with the kernel function \f$J(y,y_0)\f$.
   * \param k wave number
   * \param a scale factor
   * \param a_final scale factor at which the yukawa potential is evaluated
   */
  double sigma_J_sq (double k, double a, double a_final = 1.0);
  double sigma_J_prime_sq (double k, double a, double a_final = 1.0);

  /**
   * Returns the time derivative of \f$\sigma_J^2\f$.
   * \param k wave number
   * \param a scale factor
   * \param a_final scale factor at which the yukawa potential is evaluated
   */
/*
  double sigma_J_sq_dot (double k, double a, double a_final = 1.0);
  /**
   *Returns the derivative of \f$\sigma_J^2\f$ with respect to the amplitude of
   *the initial density fluctuation power spectrum.
   * \param k wave number
   * \param a scale factor
   * \param a_final scale factor at which the yukawa potential is evaluated
   */
  /*double dsigma_J_sq_dA (double k, double a, double a_final = 1.0);

  /**
   * Returns the averaged potential gradient
   * \f[
   *   v_{12}^\parallel := \left\langle
   *     \mathrm{i}\vec k\cdot\vec\nabla_1v_2
   *   \right\rangle
   * \f]
   * projected on the wave vector \f$\vec k\f$.
   * \param k wave number
   * \param a scale factor
   * \param a_final scale factor at which the yukawa potential is evaluated
   
  double v12_parallel (double k, double a, double a_final = 1.0);

  /**
   * Returns the averaged force
   * \f[
   *   f_{12}^\parallel := \left\langle
   *     -\mathrm{i}\vec k\cdot\vec f_{12}
   *   \right\rangle
   * \f]
   * projected on the wave vector \f$\vec k\f$.
   * \param k wave number
   * \param a scale factor
   * \param a_final scale factor at which the yukawa potential is evaluated
   
  double f12_parallel (double k, double a, double a_final = 1.0);

  /**
   * Returns the averaged interaction term
   * \f[
   *   \left\langle S_{I,12}\right\rangle = 2\int_0^t dt'\,
   *     g_{qp}(t,t')\,f_{12}^\parallel(t')
   * \f]
   * between any particle pair.
   * \param k wave number
   * \param a scale factor
   */
  double operator () (double k, double a);

  /**
   * Returns the Yukawa scale estimated directly from the interaction term.
   * \param a scale factor
   */
  double get_Yukawa_scale (double a);

  void set_factors (double f_k_in, double f_r_in, double a);

  double Dplus(double a, double k);

  double f(double a, double k);

  double E(double a);

  double m(double a);

  double g_qp_hamilton(double a, double a_prime);

  double linearP (double k, double a);

  double analyticP (double k, double a);

  double numericalP (double k, double a);
  
  double taylorP (double k, double a);
};

#endif
