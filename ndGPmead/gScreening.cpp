#include <cmath>
#include "gScreening.h"
#include <astro/cosmology/cosmologicalParameters.h>
#include "screeningParameters.h"

gScreening::gScreening (astro::cosmologyBase * cos_model_in):
cos_model (cos_model_in), k_internal (1.0)
{}

gScreening::~gScreening ()
{}
  
double gScreening::operator () (double a, double k)
{
  if (std::isnan(k) == false)
    this->k_internal = k;
  screeningParameters q (cos_model,a);
  double k_star = 10.0; /* till now, free parameter*/
  double Geff_Gn = 1.0 + q.param_B*q.param_b*pow (k_star/this->k_internal,q.param_a)*(pow(1.0 + pow(this->k_internal/k_star,q.param_a),1/q.param_b)-1.0);
  
  return Geff_Gn;
}
