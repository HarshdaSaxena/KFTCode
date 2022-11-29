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

  switch (type)
  {

   /*
    * Write power spectra to table
    */
    case (1):
    { 
      astro::parameterFile parameters ("optimised.txt");

      forceTerm SI (&D_lcdm,&cos_model_1,NULL,a_final,1);
      double k_scale = SI.get_Yukawa_scale(a_final);
      SI.set_factors
	(parameters.get("k_factor"), parameters.get("r_factor"), a_final);

      forceTerm SI2 (&D_lcdm,&cos_model_1,&g_screening,a_final,0);
      double k_factor2 = parameters.get("k_factor")*k_scale/SI2.get_Yukawa_scale(a_final);
      SI2.set_factors
	(k_factor2, parameters.get("r_factor"), a_final);

      std::ofstream outfile;
      outfile.open("temp.d");

      astro::table meadDelta ("mead.txt");

      for(int y = 1 ; y < 1000 ; y = y + 5)
        {
           for(int x = 1 ; x < 1000 ; x = x + 5)
 	{
          outfile << x/1000.0 <<","<< y/1000.0 << ","<< SI.forceTerm::derivative(y/1000.0,x/1000.0) << '\n';
        }
        }
/*
      astro::functionWriter write (output);
      write.add_header("#Column 1: ratio wavenumber y");
      write.add_header("#Column 2: scale factor x");
      write.add_header("#Column 3: Functional derivative");
      write.push_back ([&] (double y, double x) { return SI.forceTerm::derivative(y,x); });
      write (0.001, 1.0, 128, astro::LOG_SPACING);
      write (0.001,1.0,128, astro::LOG_SPACING);
*/    
      outfile.close();
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
