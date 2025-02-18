#ifndef SIDE_CREEP_H_
#define SIDE_CREEP_H_
/**
 * Function that sets up creep properties on the side of the fault
 *
 */

#include <vector>
//#include <math.h>
#include <iostream>
#include "creep_law.h"
#include "temperature_profile.h"

/* Adjust the grid activation energies and Arrhenius amplitude within
   the side creep mask to mimic the behaviour in Ben-Zion (1996) */
void sidecreepmask(
		   std::vector<double> &A, // IN/OUT: column ordered grid of Arrhenius amplitudes
		   std::vector<double> &E, // IN/OUT: column ordered grid of Activation energies on fault
                   const std::vector<double> &n, // IN: column ordered grid of fault stress exponents
                   const std::vector<double> &T, // IN: column ordered grid of fault temperatures
		   const int &nL, // Along-strike grid dimension
		   const int &nD, // Down-dip grid dimension
		   const double &cellLength, // slip cell length
		   const double &xDB ); // Creeping boundary width

void Arrhenius_amp( std::vector<double> &A,
                    const std::vector<double> &E,
                    const std::vector<double> &n,
                    const int &nL, const int &nD,
                    const double &cellLength, const double &cellHeight,
                    const double &xBD, const double &zBD,
                    const double loadingStrainRate,
                    const double sigma_eff,
					const double &Tsurface,
                    const double &Tgradient);

void Arrhenius_amp( std::vector<double> &A,
                    const std::vector<double> &E,
                    const std::vector<double> &n,
                    const std::vector<double> &T,
                    const double &A_amp,
                    const int &nL, const int &nD,
                    const double &cellLength,
                    const double &xBD);

void BZ96_boundary(std::vector<double> &A,
                   const std::vector<double> &E,
                   const std::vector<double> &n,
                   const int &nL, const int &nD,
                   const double &cellLength, const double &cellHeight,
                   const double &xBD, const double &zBD,
                   const double loadingStrainRate,
                   const double sigma_eff,
                   const double tau_ratio,
                   const double Bterm_length,
                   const double &Tsurface,
                   const double &Tgradient);


#endif /*SIDE_CREEP_H_*/
