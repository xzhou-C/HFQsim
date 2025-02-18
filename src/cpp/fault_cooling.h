#ifndef FAULT_COOLING_H_
#define FAULT_COOLING_H_

#include <math.h>
#include <vector>
#include "romberg_integration.h"
//#include "units.h"

//------------------------------------------------------------------------------
/**
    A halfspace with constant heat decay properties surrounds a fault
    with constant finite width.

    Cardwell et al (1978)
 */
class FaultCoolingModel{

public:
  FaultCoolingModel(
    const double &density, // default is 2.6 g/cm^3
    const double &diffusivity, // default is 0.01 cm^2/s
    const double &specHeat, // default is 790 J/(kg K)
    const double &faultWidth, // slip zone width, default 5 cm
    const double &cooling_distance // cooling distance, default 100 m
  ):
    _density(density),
    _diffusivity(diffusivity) , 
    _specHeat(specHeat), 
    _faultWidth(faultWidth), 
    _cooling_distance(cooling_distance) {}

public:
  /* All members are static and must be defined at in the initial part
     of the main source file.  Use static because should have only one
     decay function, and makes passing a function for integration
     easier. */

  /* Temperature generated immediately from heat rate spike*/
  double temperature(
		     const double &heatRate); // heat rate in J/

  void finite_cooling( std::vector<double> &fault_dT,
                       const double &timeStep,
                       const std::vector<double> &dq,
                       const int &nCells
                       );

  double finite_cooling_cell( const double &dT,
                              const double &dq,
                              const double &duration 
                              );  
   
   
protected:
  double _density, _diffusivity, _specHeat, _faultWidth, _cooling_distance;

};
//------------------------------------------------------------------------------
#endif /*FAULT_COOLING_H_*/
