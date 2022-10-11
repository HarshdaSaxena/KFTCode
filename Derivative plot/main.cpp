#include <iostream>
#include <fstream>
#include <string>

#include "astro/utilities/utilities.h"
#include "astro/io/clArguments.h"
#include "astro/io/functionWriter.h"
#include "astro/cosmology/cosmicStructures.h"
#include "astro/cosmology/standardGrowthFactor.h"
/*
#include "procaParameters.h"
#include "procaExpansionModel.h"
#include "procaDarkEnergyModel.h"
*/
#include "gScreening.h"

#include "forceTerm.h"

#include "screeningParameters.h"


int main (int argc, char * argv[])
{
  /**
   * Read and process command-line arguments
   */
  astro::clArguments ca (argc, argv);

  double qv    = ca.get ("-q", 2.0);
  int type     = ca.get ("-t", 1);
  double a_final = ca.get ("-a", 1.0);
  std::string model  = ca.get ("-m", "Model1");
  std::string input  = ca.get ("-i", "parameters.d");
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
  /**
   * Initialize Proca model
   */
/*
  procaParameters p_proca = setParameters (input, model);
  p_proca.qv = qv;             // Set q_V parameter
  procaExpansionModel e_proca; // Initialize Proca expansion function
*/
  /**
   * Initialize Proca dark-energy model and substitute it in Proca
   * cosmological model
   
  procaDarkEnergyModel de_proca (&e_proca, cos_model_2.getParameters ());
  cos_model_2.setDarkEnergyModel (&de_proca);
*/
  gScreening g_screening (&cos_model_1); // Initalize effective gravitational coupling
  /**
   * Initialize growth factors for both models
   */
  astro::standardGrowthFactor D_lcdm  (&cos_model_1);
  astro::standardGrowthFactor D_screening (&cos_model_1,&g_screening);
  std::ofstream outfile;
  outfile.open("temp.d");

  switch (type)
  {

   /*
    * Write tables containing the lcdm power spectrum and the exact proca spectrum.
    * Normalised at a=0.001.
    */
    case (1):
    {
      forceTerm SI (&D_lcdm,&cos_model_1,NULL,a_final,1);
      forceTerm SI2 (&D_screening,&cos_model_1,&g_screening,a_final,0);
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
