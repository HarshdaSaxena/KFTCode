#ifndef SCREENING_PARAMETERS_H
#define SCREENING_PARAMETERS_H

#include "astro/cosmology/cosmologyBase.h"

class screeningParameters
{
protected:
  astro::cosmologyBase * cos_model;
public:
  double p0, p1, p2, p3; /*Till now dont set y0 so p4,p5,p6,p7 not needed*/
  double param_a, param_b, param_B;
  double a_internal, a_min;
  screeningParameters (astro::cosmologyBase * cos_model, double a_in);
  void update_a (double a_in);
};


#endif
