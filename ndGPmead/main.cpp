#include <iostream>
#include <fstream>
#include <string>

#include "astro/utilities/utilities.h"
#include "astro/io/clArguments.h"
#include "astro/io/functionWriter.h"
#include "astro/cosmology/cosmicStructures.h"
#include "astro/cosmology/standardGrowthFactor.h"
#include "astro/utilities/nlModeller.h"
#include "gScreening.h"

#include "forceTerm.h"

#include "screeningParameters.h"


int main (int argc, char * argv[])
{
  /**
   * Read and process command-line arguments
   */
  astro::clArguments ca (argc, argv);
  int type     = ca.get ("-t", 1);
  double a_final = ca.get ("-a", 1.0);
  std::string output = ca.get ("-o", "temp.d");

  /**
   * Set standard cosmological parameters
   */
  double omega_m0 = 0.26; // dark-matter density today
  double omega_d0 = 0.7;  // dark-energy density today
  double hubble   = 0.7;  // Hubble constant in 100 km/s/Mpc today
  double omega_b0 = 0.04; // baryon density today
  /**
   * Initialize standard cosmological model
   */
  astro::cosmologyBase cos_model_1 (omega_m0, omega_d0, hubble, omega_b0);
  cos_model_1.setDarkUniverse (); // Switch off radiation density
  gScreening g_screening (&cos_model_1);

  /**
   * Initialize growth factor
   */
  astro::standardGrowthFactor D_lcdm  (&cos_model_1);
  astro::standardGrowthFactor D_screening (&cos_model_1,&g_screening);

  switch (type)
  {

   /*
    * Write power spectra to table
    */
    case (1):
    {
      forceTerm SI (&D_lcdm,&cos_model_1,NULL,a_final,1);
      forceTerm SI2 (&D_screening,&cos_model_1,&g_screening,a_final,0);

      astro::table meadDelta ("mead.txt");

      astro::functionWriter write (output);
      write.add_header("#Column 1: wavenumber k");
      write.add_header("#Column 2: non-linear GR");
      write.add_header("#Column 3: non-linear screening");
      write.add_header("#Column 4: Mead et al");
      write.add_header("#evaluated at a = "+std::to_string(a_final));
      write.push_back ([&] (double k) { return SI.analyticP(k,a_final); });
      write.push_back ([&] (double k) { return SI2.analyticP(k,a_final); });
      write.push_back ([&] (double k)
        { return meadDelta (k)*(2.0*M_PI*M_PI)/gsl_pow_3 (k); });
      write (0.01, 10.0, 128, astro::LOG_SPACING);
      break;
    }

    case (2):
    {
      forceTerm SI (&D_lcdm,&cos_model_1,NULL,a_final,1);
      forceTerm SI2 (&D_screening,&cos_model_1,&g_screening,a_final,0);

      astro::table meadDelta ("mead.txt");
      std::function<double (double)> meadEtalP = [&] (double k)
        { return meadDelta (k)*(2.0*M_PI*M_PI)/gsl_pow_3 (k); };

      const int n = 4;
      std::vector<double> k = astro::fill (astro::LOG_SPACING, n, 0.2, 5.0);
      std::vector<double> P (n);
      for (unsigned int i = 0; i < n; i++)
        P[i] = meadEtalP (k[i])/SI2.linearP (k[i], a_final);

      astro::data_model model;
      model.n = n;
      model.p = 2;
      model.x = k.data ();
      model.y = P.data ();
      model.log = false;
      model.model = [&] (double k, std::vector<double> p)
      {
        SI2.set_factors (p[0], p[1], a_final);
        return exp (SI2 (k, a_final));
      };
      model.dmodel = NULL;
      astro::nlModeller fit (model, false);
      double p[2] = {1.0, 1.0};
      fit (p);

      std::ofstream output ("optimised.txt");
      output
        << "# options containing parameters for optimal overlap"
        << std::endl
        << "k_factor  " << p[0]
        << std::endl
        << "r_factor  " << p[1]
        << std::endl;
      output.close ();

      break;
    }

    case (3):
    {
      forceTerm SI (&D_lcdm,&cos_model_1,NULL,a_final,1);
      forceTerm SI2 (&D_screening,&cos_model_1,&g_screening,a_final,0);

      astro::parameterFile parameters ("optimised.txt");
      SI2.set_factors
        (parameters.get ("k_factor"), parameters.get ("r_factor"), a_final);

      astro::table meadDelta ("mead.txt");

      astro::functionWriter write (output);
      write.push_back ([&] (double k) { return SI2.linearP (k, a_final); });
      write.push_back ([&] (double k) { return SI2.analyticP (k, a_final); });
      write.push_back ([&] (double k)
        { return meadDelta (k)*(2.0*M_PI*M_PI)/gsl_pow_3 (k); });
      write (0.01, 10.0, 128, astro::LOG_SPACING);

      break;
    }

    default:
    {
      throw std::runtime_error ("unknown option encountered, exciting!");
      break;
    }

   }


  return 0;
}
