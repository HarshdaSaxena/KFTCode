#include "astro/io/parameterFile.h"
#include "screeningParameters.h"
#include <astro/cosmology/cosmologicalParameters.h>

screeningParameters::screeningParameters (astro::cosmologyBase * cosmological_model_in, double a_in):
cos_model (cosmological_model_in), a_internal (a_in)
{
  this->a_min = 0.001;

  double r_c = 2.0/70.0;/*cosmological_model->hubbleFunction(1);*/ 
  double beta = 1.0+2.0*r_c*(cos_model->hubbleFunction(this->a_internal)/cos_model->hubbleFunction(this->a_min))*(1.0+((cos_model->dexpansionFunction_da(this->a_internal)/cos_model->dexpansionFunction_da(this->a_min))*70.0*(this->a_internal/this->a_min))/(3.0*(cos_model->hubbleFunction(this->a_internal)/cos_model->hubbleFunction(this->a_min)))) ;
 
  this->p0 = 1.0;
  this->p1 = 2.0;
  this->p2 = 1.0/(3*beta);
  this->p3 = 1.5;
  this->param_a  = this->p1 * this->p3 / (this->p1 - 1.0);
  this->param_b  = this->p1;
  this->param_B  = this->p2;
}


