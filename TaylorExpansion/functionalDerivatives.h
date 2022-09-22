#ifndef FUNCTIONAL_DERIVATIVES_H
#define FUNCTIONAL_DERIVATIVES_H


#include <astro/utilities/integrator.h>
#include <astro/utilities/functionTable.h>

#include "forceTerm.h"

class functionalDerivatives
{
protected:
  double a_min,a_max,k_min,k_max,x_min,x_max,a_final,deltaA;
  int nt;

  std::function<double(double)> deltaE;
  std::function<double(double)> deltaG;

  forceTerm * SI;

  astro::functionTable
    dDplus_da_tab,
    dDplus_dE_tab,
    dDplus_dG_tab,
    dSI_dE_tab,
    dSI_dG_tab,
    dSI_dA_tab;

public:
  /**
   * Constructor initialisting the class.
   * \param deltaE_in pointer to function describing the change in the
   * expansion function
   * \param deltaG_in pointer to function describing the change in the
   * gravitational constant
   * \param SI_in pointer to an input mean-field forceTerm
   * \param a_final_in scale factor at which the spectrum is to be evaluated
   */
  functionalDerivatives(std::function<double(double)> & deltaE_in,
  std::function<double(double)> & deltaG_in,forceTerm * SI_in,
  double a_final_in = 1.0);

  ~functionalDerivatives();

  void initFunctionTables();

  double Dplus(double a);

  double dDplus_da(double a);

  double E(double a);

  double scalefactor(double a);

  /*double dDplus_dE(double x, double a);

  double dDplus_dG(double x, double a);

  double dsigmaJsq_dG(double k, double x, double a);

  double dsigmaJsq_dE(double k, double x, double a);*/

  double dgh_dG(double x, double a, double a_prime);

  /*double dgh_dE(double x, double a, double a_prime);

  double dPlin_dE(double k,double x,double a);*/

  double dPlin_dG(double k,double x, double a);

 /* double dSI_dE(double k,double x); */

  double dSI_dG(double k,double x);

 /* double dSI_dA(double k);*/

  double dP_dG (double k, double x);

 /* double dP_dE (double k, double x);

  double A_G (double k, double x);

  double A_E (double k, double x); */
  /*
   * Returns the power spectrum for an alternative theory of gravity
   * normalized at a_ini=0.001
   */
  double P_AG(double k);
  /*
   * Returns the power spectrum for an alternative theory of gravity
   * normalized at a_final.
   */
  double operator () (double k);
  /*
   * Returns the lcdm reference power spectrum.
   */
  double P_GR(double k);

};

#endif
