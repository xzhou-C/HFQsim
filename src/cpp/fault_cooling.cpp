/**
    Source file for the member functions of FaultCoolingModel
*/
#include "fault_cooling.h"

//------------------------------------------------------------------------------
/* Imediate temperature from heat rate over a certain area */         //there is another area in the volume term
double FaultCoolingModel::temperature( const double &heatRate)
{
  return 0.5*heatRate/( _faultWidth * _density * _specHeat ); 
}
                                  
void FaultCoolingModel::finite_cooling( std::vector<double> &fault_dT,
                                        const double &timeStep,
                                        const std::vector<double> &dq,
                                        const int &nCells)
{
  #pragma omp parallel
  {
    #pragma omp for
    for( unsigned int k=0; k<nCells; k++){
      fault_dT[k] = finite_cooling_cell( fault_dT[k], dq[k], timeStep);
    }
  }
  return;
}

double FaultCoolingModel::finite_cooling_cell( const double &dT,
                                               const double &dq,
                                               const double &duration ) 
{
  return (dT + temperature(dq) ) * exp(- _diffusivity / (_faultWidth * _cooling_distance) * duration );
}
