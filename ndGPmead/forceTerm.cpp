#include "forceTerm.h"
#include "screeningParameters.h"
#include <astro/utilities/nlModeller.h>
#include <KFT/particleDynamics/zeldovichParticleDynamics.h>
#include <KFT/particleDynamics/comovingParticleDynamics.h>
#include <astro/cosmology/gNewton.h>

forceTerm::forceTerm
  (astro::abstractGrowthFactor * growth_model_in,
   astro::cosmologyBase * cosmological_model_in,
   astro::gEffective * g_effective_in, double a_norm_in, int proptype):
a_min (0.001), a_max (1.0), k_min (0.1), k_max (10.0), k_scale (0.0), k_star (1.0),
a_save (-1.0), k_save (-1.0), f_k_save (-1.0), f_r_save (-1.0), n_table (128),
g_effective (g_effective_in), proptype (proptype)
{
  if (cosmological_model_in != NULL)
  {
      cos_model_internal = false;
      cosmological_model = cosmological_model_in;
  }
  else
  {
      cos_model_internal = true;
      cosmological_model = new astro::cosmologyBase;
      cosmological_model->setDarkUniverse ();
  }

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

  propagator = new KFT::zeldovichParticleDynamics
    (cosmological_model, a_min, growth_model_in);

  propagator_h = new KFT::comovingParticleDynamics
    (cosmological_model, a_min, growth_model_in);

  tophat_filter = new astro::tophatFilter;
  gauss_filter  = new astro::gaussFilter;

  linear_spectrum = new astro::bardeenPowerSpectrum
    (cosmological_model, 0.8, tophat_filter);
  nonlinear_spectrum = new astro::smithEtalPowerSpectrum
    (cosmological_model, 0.8, tophat_filter);

  cosmic_structures = new astro::cosmicStructures (cosmological_model);
  if (growth_model_in != NULL)
      cosmic_structures->setGrowthFactor (growth_model_in);

  proptype = proptype;
  a_norm = a_norm_in;

  D_plus_min = cosmic_structures->Dplus (a_min);
  D_plus_max = cosmic_structures->Dplus (a_max);
  D_plus_norm = cosmic_structures->Dplus (a_norm);
  t_max = D_plus_max/D_plus_min-1.0;
  amplitude = D_plus_min*D_plus_min/D_plus_norm/D_plus_norm;
  sigma_1 = sqrt (amplitude*linear_spectrum->sigma2 (0.0, a_norm, 1));

  viscosity = gsl_pow_2 (linear_spectrum->nonlinearRadius (1.0/D_plus_max));

  init_tau ();
  init_tables ();
  std::cout<<tau<<std::endl; 

}

forceTerm::~forceTerm()
{
  if (g_eff_internal)
    delete g_effective;
  if (cos_model_internal)
    delete cosmological_model;

  delete propagator;
  delete propagator_h;
  delete tophat_filter;
  delete gauss_filter;
  delete linear_spectrum;
  delete nonlinear_spectrum;
  delete cosmic_structures;
}

void forceTerm::init_tau ()
{
  linear_spectrum->setFilter (gauss_filter);

  std::function<double (double)> dx = [&] (double t)
  {
    astro::integrator integrate
      ([&] (double tau)
        {
          return
            sqrt (amplitude*
            linear_spectrum->sigma2 (sqrt (2.0*viscosity*tau), a_norm, 1))/
            sigma_1;
        });
    return integrate (0.0, t);
  };
  std::vector<double>
    t = astro::fill (astro::LOG_SPACING, n_table, 0.1*t_max, t_max);
  std::vector<double> x (n_table);
  for (unsigned int i = 0; i < n_table; i++)
    x[i] = dx (t[i]);

  linear_spectrum->setFilter (tophat_filter);

  std::function<double (double, std::vector<double>)>
    mf = [&] (double t, std::vector<double> p)
    {
      tau = p[0];
      return displacement (t);
    };
  std::function<std::vector<double> (double, std::vector<double>)>
    dm = [mf] (double x, std::vector<double> p)
    {
      std::vector<double> df (1);
      df[0] = gsl_pow_2 (mf (x, p))/(2.0*p[0]*sqrt (x*p[0]));
      return df;
    };
  astro::data_model model;
  model.n = n_table;
  model.p = 1;
  model.x = t.data ();
  model.y = x.data ();
  model.log = false;
  model.model = mf;
  model.dmodel = dm;

  astro::nlModeller nl_modeller (model);
  double p[1] = {25.0};
  nl_modeller (p);
  tau = p[0];
}

void forceTerm::init_tables (double f_k, double f_r)
{
  bool changed = (fabs (k_scale-k_save) > 1.0e-4);
  if ((! changed) && (k_save > 0.0))
    return;

  k_save = k_scale;

  sigma_J_sq_table.fill
    ([&] (double k, double a)
    {
      double r = sigma_1*f_r*displacement (std::max(0.0,propagator->g_qp (a)));
      astro::integrator integrate ([&] (double y)
      { return
          y*y*amplitude*linear_spectrum->operator () (y,a_norm)/
          (1.0+y*y*r*r)*
          J (y/k, f_k*k_scale/k)/gsl_pow_2 (2.0*M_PI); });
      return integrate (0.0);
    }, astro::LOG_SPACING, astro::LOG_SPACING, n_table, n_table,
    k_min, k_max, a_min, a_max);

  sigma_J_prime_sq_table.fill
    ([&] (double k, double a)
    {
      double r = sigma_1*f_r*displacement (std::max(0.0,propagator->g_qp (a))); /* r = sigma_v^2 lambda^2*/
      astro::integrator integrate ([&] (double y)
      { return
          y*y*amplitude*linear_spectrum->operator () (y,a_norm)/
          (1.0+y*y*r*r)*
          J_prime (y/k, k_scale/k, k_star/k)/gsl_pow_2 (2.0*M_PI); },1.0e-2); /* integrate y^2 linear spectrum == P_i by damping, y_0/k = k_scale/k, this is the changed sigma_prime*/
      /*std::cout<<integrate (0.0);*/
      return integrate (0.0);
    }, astro::LOG_SPACING, astro::LOG_SPACING, n_table, n_table,
    k_min, k_max, a_min, a_max);
/*
  sigma_J_sq_dot_table.fill
    ([&] (double k, double a)
    {
      double r = sigma_1*f_r*displacement (std::max(0.0,propagator->g_qp (a)));
      astro::integrator integrate ([&] (double y)
      { return
          -y*y*y*y*r*r/pow((1.0+y*y*r*r),2.0)
          *(2.0/propagator->g_qp (a)-1.0/
          (sqrt(tau*std::max(0.0,propagator->g_qp (a)))+propagator->g_qp (a)))
          *amplitude*linear_spectrum->operator () (y,a_norm)*
          J (y/k, f_k*k_scale/k)/gsl_pow_2 (2.0*M_PI); });
      return integrate (0.0);
    }, astro::LOG_SPACING, astro::LOG_SPACING, n_table, n_table,
    k_min, k_max, a_min, a_max);

  dsigma_J_sq_dA_table.fill
    ([&] (double k, double a)
    {
      double r = sigma_1*f_r*displacement (std::max(0.0,propagator->g_qp (a)));
      astro::integrator integrate ([&] (double y)
      { return
          y*y*amplitude*linear_spectrum->operator () (y,a_norm)/
          (1.0+y*y*r*r)*(1.0-y*y*r*r/(1.0+y*y*r*r))*
          J (y/k, f_k*k_scale/k)/gsl_pow_2 (2.0*M_PI); });
      return integrate (0.0);
    }, astro::LOG_SPACING, astro::LOG_SPACING, n_table, n_table,
    k_min, k_max, a_min, a_max);

    if(proptype==1)
      init_tables_zeldovich(f_k, f_r);
*/

}

/*
void forceTerm::init_tables_zeldovich (double f_k, double f_r)
{
  v12_parallel_table.fill
    ([&] (double k, double a)
    {
      return
        g_effective->operator()(a)*propagator->g_dot (a)/propagator->g (a)*
        cosmic_structures->Dplus (a)/D_plus_min*sigma_J_sq_table (k, a);

    }, astro::LOG_SPACING, astro::LOG_SPACING, n_table, n_table,
    k_min, k_max, a_min, a_max);

  f12_parallel_table.fill
    ([&] (double k, double a)
    {
      astro::integrator integrate ([&] (double x)
        { return propagator->g (x)*v12_parallel_table (k, x)*
          propagator->jacobian (x)/x; });
      return
        v12_parallel_table (k, a)-
        g_effective->operator()(a)*propagator->g_dot (a)/
        gsl_pow_2 (propagator->g (a))*integrate (a_min, a);

    }, astro::LOG_SPACING, astro::LOG_SPACING, n_table, n_table,
    k_min, k_max, a_min, a_max);

}
*/

double forceTerm::displacement (double t)
{
  return t/(1.0+sqrt (t/tau)); }

//include asymtotic behaviour for large arguments to prefent oscilations
//in sigma_J_sq_dot
double forceTerm::J (double y, double y0)
{
  return (y<=100)? (1.0+0.25*(1.0-y*y-y0*y0)/y*
  log ((y0*y0+gsl_pow_2 (1.0+y))/(y0*y0+gsl_pow_2 (1.0-y))))
  : (2.0/3.0/gsl_pow_2(y));
}

double forceTerm::J_prime (double y, double y0, double y_star) /*J_prime definition*/
{ screeningParameters q (cosmological_model, 0.001);
  astro::integrator integrate_J_prime ([&] (double x)
      { return
          (1-y*x)/(1+y0*y0 + y*y -2*x*y)/pow(1+y*y-2*x*y,q.param_a/2) *(pow(1+ pow(1+y*y-2*y*x,q.param_a/2)/pow(y_star,q.param_a),1/q.param_b)-1); },1.0e-2);
  /*std::cout<<integrate_J_prime (-1.0,1.0);*/
  return integrate_J_prime (-1.0,1.0);
}


double forceTerm::sigma_J_sq (double k, double a, double a_final)
{
    if (fabs (a_final-a_save) > 1.0e-3)
    {
      k_scale = get_Yukawa_scale (a_final);
      init_tables();
    }
    a_save = a_final;

  return sigma_J_sq_table (k, a); }

double forceTerm::sigma_J_prime_sq (double k, double a, double a_final)
{
    if (fabs (a_final-a_save) > 1.0e-3)
    {
      k_scale = get_Yukawa_scale (1.0);
      init_tables();
    }
    a_save = a_final;

  return sigma_J_prime_sq_table (k, a); }
/*

double forceTerm::sigma_J_sq_dot (double k, double a, double a_final)
{
    if (fabs (a_final-a_save) > 1.0e-3)
    {
      k_scale = get_Yukawa_scale (a_final);
      init_tables();
    }
    a_save = a_final;

  return sigma_J_sq_dot_table (k, a); }

double forceTerm::dsigma_J_sq_dA (double k, double a, double a_final)
{
    if (fabs (a_final-a_save) > 1.0e-3)
    {
      k_scale = get_Yukawa_scale (a_final);
      init_tables();
    }
    a_save = a_final;

  return dsigma_J_sq_dA_table (k, a); }

double forceTerm::v12_parallel (double k, double a, double a_final)
{
    if (fabs (a_final-a_save) > 1.0e-3)
    {
      k_scale = get_Yukawa_scale (a_final);
      init_tables();
    }
    a_save = a_final;

  return v12_parallel_table (k, a); }

double forceTerm::f12_parallel(double k, double a, double a_final)
{
    if (fabs (a_final-a_save) > 1.0e-3)
    {
      k_scale = get_Yukawa_scale (a_final);
      init_tables();
    }
    a_save = a_final;

   return f12_parallel_table (k, a); }
*/

double forceTerm::interaction_term_hamilton (double k, double a) /*3*int g_H *G_eff* 1/(a^2 D_plus' E)* D_plus^2*sigma_J^2 D_plus' gets cancelled from t->a variable change in intg*/
{ screeningParameters q (cosmological_model,a);
  astro::integrator integrate
    ([&] (double x)
    { //astro::cosmologicalParameters * p = cos_model->getParameters ();
      //screeningParameters * q = reinterpret_cast<screeningParameters*> (p->par);
      return
        1./a_min*propagator_h->g_qp(a,x)*
        pow(Dplus(x,k),2.)*(sigma_J_sq_table(k,x)+ q.param_B*q.param_b*sigma_J_prime_sq_table(k,x)*pow(k_star/k,q.param_a))/pow(x/a_min,2.)/E(x); },1.0e-2);
  /*std::cout<<integrate (a_min, a);*/
  return 3.0*integrate (a_min, a);
}

double forceTerm::interaction_term_normal (double k, double a) /*3*int g_H *G_eff* 1/(a^2 D_plus' E)* D_plus^2*sigma_J^2 D_plus' gets cancelled from t->a variable change in intg*/
{ astro::integrator integrate
    ([&] (double x)
    { //astro::cosmologicalParameters * p = cos_model->getParameters ();
      //screeningParameters * q = reinterpret_cast<screeningParameters*> (p->par);
      return
        1./a_min*propagator_h->g_qp(a,x)*
        pow(Dplus(x,k),2.)*(sigma_J_sq_table(k,x))/pow(x/a_min,2.)/E(x); },1.0e-2);
  /*std::cout<<integrate (a_min, a);*/
  return 3.0*integrate (a_min, a);
}

double forceTerm::operator () (double k, double a)
{
    if (fabs (a-a_save) > 1.0e-3)
    {
      k_scale = get_Yukawa_scale (1.0);
      init_tables();
    }
    a_save = a;

  switch(proptype){
    case (0):{
      return interaction_term_hamilton (k, a);
      break;}
    case(1):{
      return interaction_term_normal (k, a);
      break;}
    default:{
      throw std::runtime_error ("Error: Unknown propagator type");
      break;}
  }

}

double forceTerm::get_Yukawa_scale (double a)
{
  std::function<double (double, std::vector<double>)>
    model = [] (double x, std::vector<double> p)
    { return p[0]*x*x/(p[1]*p[1]+x*x); };

  std::function<std::vector<double> (double, std::vector<double>)>
    dmodel = [&] (double x, std::vector<double> p)
    {
      std::vector<double> df;
      df.push_back (model (x, p)/p[0]);
      df.push_back (-2.0*p[1]*model (x, p)/(p[1]*p[1]+x*x));
      return df;
    };

  astro::data_model data;
  data.n = n_table;
  data.p = 2;
  data.x = new double[data.n];
  data.y = new double[data.n];
  data.log = false;
  data.model  = model;
  data.dmodel = dmodel;

  for (unsigned int i = 0; i < data.n; i++)
  {
    data.x[i] = astro::x_linear (i, data.n, 0.05, 5.0); /*f_v*/

    switch(proptype){
      case (0):{
      data.y[i] = 1.0-exp (-0.5*interaction_term_hamilton (data.x[i], a));
      break;}
      case(1):{
      data.y[i] = 1.0-exp (-0.5*interaction_term_normal (data.x[i], a));
      break;}
      default:{
        throw std::runtime_error ("Error: Unknown propagator type");
        break;}

      }
  }

  astro::nlModeller nl_model (data);
  double p[2] = {1.0, 1.0};
  nl_model (p);
  return p[1];
}

void forceTerm::set_factors (double f_k, double f_r, double a)
{
  k_scale = 0.0;
  init_tables ();
  k_scale = get_Yukawa_scale (a);
  init_tables (f_k, f_r);

  a_save = a;
  f_k_save = f_k;
  f_r_save = f_r;
}

double forceTerm::Dplus(double a, double k)
{ return cosmic_structures->Dplus(a,k)/D_plus_min;}

double forceTerm::f(double a, double k)
{ return cosmic_structures->dlnGrowthFactor_dlna(a)/
  cosmic_structures->dlnGrowthFactor_dlna(a_min);}

double forceTerm::E(double a)
{ return cosmological_model->expansionFunction(a)/
  cosmological_model->expansionFunction(a_min);}

/*double forceTerm::g_qp_zeldovich (double a, double a_prime)
{ return propagator->g_qp(a);}*/

double forceTerm::g_qp_hamilton (double a, double a_prime)
{ return propagator_h->g_qp(a,a_prime); }

double forceTerm::m(double a)
{ return propagator->g(a);}

/*double forceTerm::m_dot(double a)
{ return propagator->g_dot(a);}*/

double forceTerm::linearP (double k, double a)
{ return linear_spectrum->operator () (k, a); }

double forceTerm::analyticP (double k, double a)
{ return Dplus(a,k)*Dplus(a,k)*amplitude*linear_spectrum->operator () (k, a_norm)*
  exp (operator () (k, a));}

double forceTerm::numericalP (double k, double a)
{ return nonlinear_spectrum->operator () (k, a); }
