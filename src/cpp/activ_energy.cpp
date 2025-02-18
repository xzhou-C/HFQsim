#include "activ_energy.h"

void FaultActivationEnergy::setAE_linear( std::vector<double> &E,
                   const double &E0,                // surface activation energy
                   const double &Emax,              // maximum activation energy
                   const double &zE,                // transition depth, E = Emax below zE
                   const int &nL, const int &nD,    // fault dimension
                   const double &cellHeight)
{
    if( E.size() != nL * nD){
        std::cerr << " activation energy ERROR: size of E must be nL*nD!"<< std::endl;
        std::terminate();
    }
    // E.clear();

    int jmax = floor(zE/cellHeight-0.5);
    for( int i=0; i<nL; i++){
        for( int j=0; j<nD; j++){
            int fluct = rand() % 10 + 1;
            if( j <= jmax){
                E[i*nD+j] = E0 + (Emax - E0)*(double)j/(double)jmax;
            }
            else{
                E[i*nD+j] = Emax;
           }
        }
    }
    return;
}

void FaultActivationEnergy::setAE_exp( std::vector<double> &E,
                const double &E0,           // surface activation energy
                const double &Emax,         // maximum activation energy
                const double &zE,           // transition depth, E = Emax below zE
                // const double &dE,           // decay rate of E
                const int &nL,const int &nD,
                const double &cellHeight)
{
    if( E.size() != nL * nD){
        std::cerr << " activation energy ERROR: size of E must be nL*nD!"<< std::endl;
        std::terminate();
    }

    // int jmax = floor(zE/cellHeight-0.5);
    for( int j=0; j<nD; j++){
        double depth = cellHeight * (j + 0.5);
        for( int i=0; i<nL; i++){
            if( depth <= zE){
                E[i*nD+j] = Emax - (Emax - E0)*exp(-10*depth/zE);
            }
            else{
                E[i*nD+j] = Emax;
            }
        }
    }
    return;
}

/* add asperity */
// void AEasp_square( std::vector<double> &E,
//                 const int &TL, const int &BL,   // index of top and bottom left conor
//                 const int &TR, const int &BR)   // index of top and bottom right conor

void FaultActivationEnergy::AEasp_square( std::vector<double> &E,
                                          const int &leftX,
                                          const int &rightX,
                                          const int &topZ,
                                          const int &bottomZ,
                                          const double &ratio,
                                          const int &nL,
                                          const int &nD)
{
    if( E.size() != nL * nD){
        std::cerr << " activation energy ERROR: size of E must be nL*nD!"<< std::endl;
        std::terminate();
    }
    for( int i=leftX; i<rightX+1; i++){
        for( int j=topZ; j<bottomZ+1; j++){
            double fluct = rand() % 10 + 1;
            E[i*nD+j] = ((fluct - 5)/100 + ratio) * E[i*nD+j];
        }
    }
    return;
}

void FaultActivationEnergy::AEasp_diamond( std::vector<double> &E,
                    const int &z,
                    const int &xmin, const int &xmax,
                    const double &ratio,
                    const int &nL, const int &nD)
{
    if( E.size() != nL * nD){
        std::cerr << " activation energy ERROR: size of E must be nL*nD!"<< std::endl;
        std::terminate();
    }
    for( int i=xmin; i<xmax+1; i++){
        int dz = std::min(i-xmin,xmax-i);
        for( int j=z-dz; j<z+1+dz; j++){
            E[i*nD+j] = ratio * E[i*nD+j];
        }
    }
    return;
}

// void AEasp_circle

/* apply healing mechanism */
/*
double FaultActivationEnergy::computeAE( const double &E0,                // cellular initial activation energy
                                         const double &tNow,              // current model time
                                         const double &t_seismicSlip,     // time of last EQ
                                         const double &tau_max,           // maximum healing time
                                         const double &tau_min,           // minimum healing time
                                         const double &scale_factor)      // scaling factor of the relaxation function
{
    return E0 * ( 1.0 + scale_factor * integRelax(tNow-t_seismicSlip,tau_max,tau_min) );
}

double FaultActivationEnergy::integRelax( const double &tHeal,        // tNow - t_seismicSlip
                                          const double &tau_max,      // maximum healing time
                                          const double &tau_min)      // minimum healing time
{
    if( tHeal < 0 || tHeal >= 5*tau_max){
        return 0.0;
    }
    else{
        return R.integrate( &relaxation,tau_min,tau_max);
    }
}

double FaultActivationEnergy::relaxation( double healingTime)
{
    return (1/healingTime)*exp(-timeSinceDamage/healingTime);
}
*/
/* healing function from Zöller et al., 2004 */

void FaultActivationEnergy::logHeal( std::vector<double> &E,            // activation energy on the fault
                                     const std::vector<double> &E0,     // activation energy after EQ
                                     const std::vector<double> &Emax,   // initial activation energy
                                     // const std::vector<double> &ratio,  // scale factor
                                     const std::vector<double> &tEQ,    // time of last seismic slip at each cell
                                     const std::vector<double> &tHeal,  // healing time, depth (temperature/stress) dependent
                                     const double &tNow,                // current time
                                     const int &nL, const int &nD)
{
    if( E.size() != nL * nD){
        std::cerr << " activation energy ERROR: size of E must be nL*nD!"<< std::endl;
        std::terminate();
    }
    double ratio = 1 / log(2);
    #pragma omp parallel
    {
    #pragma omp for
    for( unsigned int i=0; i<nL; i++){
        for (unsigned int j=0; j<nD; j++){
            int iCell = i*nD+j;
            if( tEQ[iCell] <= 0.0 || tNow-tEQ[iCell] >= tHeal[iCell]){
                E[iCell] = Emax[iCell];
            }
            else{
                E[iCell] = E0[iCell] + ratio * (Emax[iCell] - E0[iCell]) * log(1+(tNow-tEQ[iCell])/tHeal[iCell]);
            }
        }
    }
    }
    return;
}
