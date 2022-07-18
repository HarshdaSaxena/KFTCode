#include <cmath>
#include "gScreening.h"
#include <astro/cosmology/cosmologicalParameters.h>

gScreening::gScreening (astro::cosmologyBase * cos_model_in):
cos_model (cos_model_in), k_internal (1.0)
{}

gScreening::~gScreening ()
{}
  
double gScreening::operator () (double a, double k)
{
  if (std::isnan(k) == false)
    this->k_internal = k;
  double r_c = 2.0/70.0; /* Using the scale of H_0rc = 2 */
  double beta = 1.0+2.0*r_c*(cos_model->hubbleFunction(a)/cos_model->hubbleFunction(0.001))*(1.0+((cos_model->dexpansionFunction_da(a)/cos_model->dexpansionFunction_da(0.001))*70.0*(a/0.001))/(3.0*(cos_model->hubbleFunction(a)/cos_model->hubbleFunction(0.001))));
  double k_star = 10.0; /* till now, free parameter*/
  double Geff_Gn = 1.0 + (2.0/(3.0*beta))*pow (k_star/this->k_internal,3.0)*(pow(1.0 + pow(this->k_internal/k_star,3.0),0.5)-1.0);
  
  return Geff_Gn;
}
