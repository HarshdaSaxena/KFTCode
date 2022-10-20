#ifndef SCREENING_G_EFFECTIVE_H
#define SCREENING_G_EFFECTIVE_H

#include "astro/cosmology/cosmologyBase.h"
#include "astro/cosmology/gEffective.h"
#include <limits>
#include "screeningParameters.h"

class gScreening: public astro::gEffective
{
protected:
  astro::cosmologyBase * cos_model;
public:
  double k_internal;
  gScreening (astro::cosmologyBase * cos_model_in);

  ~gScreening ();
  
  double operator () (double a, double k=std::numeric_limits<double>::quiet_NaN());
};

#endif
