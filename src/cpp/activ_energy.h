#ifndef ACTIV_ENERGY_H_
#define ACTIV_ENERGY_H_

#include <vector>
#include <iostream>
#include <math.h> 
#include <stdlib.h> // rand()
#include "romberg_integration.h"

class FaultActivationEnergy{

public:
    static double timeSinceDamage; // used for integral, need to be updated for each cell at each step

    void setAE_linear( std::vector<double> &E,
                       const double &E0,                // surface activation energy
                       const double &Emax,              // maximum activation energy
                       const double &zE,                // transition depth, E = Emax below zE
                       const int &nL, const int &nD,    // fault dimension
                       const double &cellHeight);

    void setAE_exp( std::vector<double> &E,
                    const double &E0,           // surface activation energy
                    const double &Emax,         // maximum activation energy
                    const double &zE,           // transition depth, E = Emax below zE
                    // const double &dE,           // decay rate of E
                    const int &nL,const int &nD,
                    const double &cellHeight);

    void AEasp_square( std::vector<double> &E,  // activation energy on the fault
                       const int &leftX,        // left boundary of the box asperity
                       const int &rightX,       // right side
                       const int &topZ,         // top side
                       const int &bottomZ,      // bottom side
                       const double &ratio,     // new_E = ratio * old_E
                       const int &nL,
                       const int &nD);

    void AEasp_diamond( std::vector<double> &E,
                        const int &z,                           // vertical index of the middle row
                        const int &xmin, const int &xmax,       // left and right edges of the asperity
                        const double &ratio,                    // new_E = ratio * old_E
                        const int &nL, const int &nD);          // fault dimension

/*    double computeAE( const double &E0,                // cellular initial activation energy
                      const double &tNow,              // current model time
                      const double &t_seismicSlip,     // time of last EQ
                      const double &tau_max,           // maximum healing time
                      const double &tau_min,           // minimum healing time
                      const double &scale_factor);     // scaling factor of the relaxation function

    double integRelax( const double &tHeal,        // tNow - t_seismicSlip
                       const double &tau_max,      // maximum healing time
                       const double &tau_min);     // minimum healing time
*/
    void logHeal( std::vector<double> &E,            // activation energy on the fault
                  const std::vector<double> &E0,     // activation energy after EQ
                  const std::vector<double> &Emax,   // initial activation energy 
                  // const std::vector<double> &ratio,  // scale factor
                  const std::vector<double> &tEQ,    // time of last seismic slip at each cell
                  const std::vector<double> &tHeal,  // healing time, depth (temperature/stress) dependent
                  const double &tNow,                // current time
                  const int &nL, const int &nD);

 private:
     Romberg R;

   //  static double relaxation( double healingTime);


 }; // end of activ_energy.h
#endif
