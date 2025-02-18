#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#include <vector>
#include <string.h>
#include <math.h>
#include "units.h"


//----------------------------------------
/**  Set parameters used in main function
*/
namespace In
{
  /* Fault properties */
  const double faultWidth( 1.0*Units::km ), // width of the fault zone for strain calc, 10 cm is the initial value.
    rigidity( 30*Units::GPa ); // rigidity of the halfspace, value for granite is 24 GPa. 

  /* Loading and creep parameters */
  const double plateVelocity( 35*Units::mm/Units::year ), // plate velocity
    loadingStrainRate( plateVelocity/faultWidth ) , // strain rate of fault loading
    // du_max( 100.0*Units::cm ), // maximum creep distance
    arrhAmplitude( 6.31E-24 / ( pow(Units::Pa, 3) * Units::second) ), // Amplitude in the Arrhenius relation  
    // srufaceEdamaged( 50*Units::kJ/Units::mole), // not using#################   
    AEasp_ratio( 0.1), // scaling factor of the asperity activation energy    
    stressExponent( 3.0 ), // exponent to stress in arrhenius equation
    xBD( 7.5*Units::km ), // creeping boundary width
    zBD( 10*Units::km), // BD boundary depth 
    zShallow( 3*Units::km), //shallow depth 
    tau_ratio( 4.0 ),
    Bterm_length( 7.5*Units::km); // parameters used in BZ96

  /* Background temperature profile */
  const double Tsurface( 278.15*Units::K ),  // Temperature at surface, 273.15 is the original value
    Tgradient( 20.0 * Units::K/Units::km ); // Change in temperature with depth

  /*  Thermal parameters */
  const double density = 2.6*Units::g/(Units::cm*Units::cm*Units::cm), // density of crustal rock
    diffusivity = 1.0E-2*Units::cm*Units::cm/Units::second, // diffusivity of granite
    specHeat = 790*Units::J/(Units::kg*Units::K); // specific heat capacity for granite
   
  /* Strength and initial stress parameters */
  const double
    fs(0.75), // coefficient of friction
    dsigmadz( 18.0 * Units::MPa/Units::km ), // effective normal stress gradient
    dsigmadz_dry( 28.0 * Units::MPa/Units::km ), // contact zone, zero pore pressure 
    dynStrengthDrop( 2.0*Units::MPa ); // tau_static - tau_dynamic

  /* EQ slip velocity */
  const double slipVelocity(1.0 * Units::m/Units::second);   

  /* Algorithm time constraints */
  const double maxTimeStep = 3.0*Units::day;
  const double minTimeStep = 0.001*Units::second;        
  const int nTimeMax = 1e9; // maximum number of time steps
  const double recordTime = 0.1*Units::year; // time to start applying heating and healing
  // const double coolTime = 0.5*Units::second; // quench temperature with fluid cooling, reset dT to zero

  /* Output file names */
  const std::string ofilesuffix_initStress = "_initstress.txt";
  const std::string ofilesuffix_finalStress = "_finalstress.txt";
  const std::string ofilesuffix_creepSlip = "_creepslip.txt";
  const std::string ofilesuffix_eqSlip = "_eqkslip.txt";
  const std::string ofilesuffix_slipdef = "_finalslipdef.txt";
  const std::string ofilesuffix_faultBkgdT = "_faultBkgdT.txt";
  const std::string ofilesuffix_faultSlip = "_faultSlip.txt";
  const std::string ofilesuffix_hist_dT = "_hist_dT.txt";
  const std::string ofilesuffix_slip_type = "_slip_type.txt";
  const std::string ofilesuffix_timeStep = "_timeStep.txt";
  const std::string ofilesuffix_hist_tau = "_hist_tau.txt";
  const std::string ofilesuffix_hist_vc = "_hist_vc.txt";
  //const std::string ofilesuffix_stressavg = "_stressavg.txt";
  //const std::string ofilesuffix_stressdev = "_stressdev.txt";
  const std::string ofilesuffix_activE = "_activE.txt";
  const std::string ofilesuffix_Arrheamp = "_Arrheamp.txt";
  const std::string ofilesuffix_StrengthChange = "_SCP.txt";
  const std::string ofilesuffix_hist_actE = "_hist_actE.txt";
  const std::string ofilesuffix_tau_heal = "_tau_heal.txt";
  const std::string ofilesuffix_slipdefHist = "_hist_slipdef.txt";
  const std::string ofilesuffix_hist_taus = "_hist_taus.txt";

  const std::string ofilesuffix_MHPStress = "_MHPstress.txt";
  const std::string ofilesuffix_MHPTemp = "_MHPtemp.txt";
  const std::string ofilesuffix_MHPVcreep = "_MHPvcreep.txt";
}

#endif /* PARAMETERS_H_ */
