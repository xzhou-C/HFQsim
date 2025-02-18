#include "side_creep.h"

//------------------------------------------------------------------------------
/**
   Adjust Arrhenius amplitude and activation energy within a distance xBD from
   the side of the fault so that...
   1. there is no temperature dependent creep for x<xBD and x>(L-xBD)
   2. the strength profile as a function of x reflects the strength profile as a
      function of (W-z) at x=xBD, where W is the fault width, and z is the depth.
*/
void sidecreepmask( std::vector<double> &A,
                    std::vector<double> &E,
                    const std::vector<double> &n,
                    const std::vector<double> &T,
		                const int &nL, const int &nD,
                    const double &cellLength,
                    const double &xBD )
{
  /* Check dimensions of input arguments */
  if( A.size() != nL*nD || E.size() != nL*nD ){
    std::cerr << "ERROR: size of A and E must be nL*nD" << std::endl;
    std::terminate();
  }

  int iSide = 0;
  while( ((double)iSide + 0.5)*cellLength < xBD ) iSide++;

  /* Loop through along strike grid coordinates each grid point */
  for( int i=0; i<iSide; i++ )
    {
      //double sideDist = ((double)i + 0.5)*cellLength;
      int iRev = nL-i-1;

      /* Get the reference cell to mimic */
      int kRef = iSide*nD + nD-i-1;
      CreepLaw C1left( A[kRef], E[kRef], n[kRef] );
      double A1left = C1left.temperatureTerm( T[kRef] );

      /* Repeat for the right hand side */
      kRef = (nL-1-iSide)*nD + nD-i-1;
      CreepLaw C1right( A[kRef], E[kRef], n[kRef] );
      double A1right = C1right.temperatureTerm( T[kRef] );

      /* Loop through down-dip grid coordinates */
      for( int j=0; j<nD; j++ )
        {
          /* Set the alternative A to mimic what is there */
          CreepLaw C2left( A[i*nD+j], E[i*nD+j], n[kRef] );
          double A2left = C2left.temperatureTerm( T[i*nD+j] );

          CreepLaw C2right( A[iRev*nD+j], E[iRev*nD+j], n[kRef] );
          double A2right = C2right.temperatureTerm( T[iRev*nD+j] );

          /* Turn off temperature dependence by setting E=0 */
          E[i*nD+j] = 0.0;
          E[iRev*nD+j] = 0.0;

          /* Choose whichever is weaker */
          A[i*nD+j] = std::max( A1left, A2left );
          A[iRev*nD+j] = std::max( A1right, A2right );

        } // END Loop through depth

    } // END Loop through bottom

  return;
}

void Arrhenius_amp( std::vector<double> &A,
                    const std::vector<double> &E,
                    const std::vector<double> &n,
                    const int &nL, const int &nD,
                    const double &cellLength, const double &cellHeight,
                    const double &xBD, const double &zBD,
                    const double loadingStrainRate,
                    const double sigma_eff,
                    const double &Tsurface,
                    const double &Tgradient)
{
  if( A.size() != nL*nD || E.size() != nL*nD ){
    std::cerr << "side creep ERROR: size of A and E must be nL*nD" << std::endl;
    std::terminate();
  }

  BackgroundTemperatureProfile T( Tsurface, Tgradient);

  /* brittle section */
  for( int j=0; ((double)j+0.5)*cellHeight < zBD; j++)
    {
      int iBrittle = floor(xBD/cellLength-0.5);
      double depth = ((double)j+0.5)*cellHeight;
      for( int i=iBrittle; i<nL-iBrittle-1; i++)
      {
        A[i*nD+j] = loadingStrainRate/(pow(sigma_eff,3)*exp(-E[i*nD+j]/(PhysicalConstants::R_g * T(depth) )));
      }
    }
  
  return;
} //End of Arrhenius_amp

void Arrhenius_amp( std::vector<double> &A,
                    const std::vector<double> &E,
                    const std::vector<double> &n,
                    const std::vector<double> &T,
                    const double &A_amp,
                    const int &nL, const int &nD,
                    const double &cellLength,
                    const double &xBD)
{
  if( A.size() != nL*nD || E.size() != nL*nD ){
    std::cerr << "side creep ERROR: size of A and E must be nL*nD" << std::endl;
    std::terminate();
  }

  int iSide = 0;
  while( ((double)iSide + 0.5)*cellLength < xBD ) iSide++;

  for( int i=0; i<iSide; i++ )
    {
      //double sideDist = ((double)i + 0.5)*cellLength;
      int iRev = nL-i-1;

      /* Get the reference cell to mimic */
      int kRef = iSide*nD + nD-i-1;
      CreepLaw C1left( A[kRef], E[kRef], n[kRef] );
      double A1left = C1left.temperatureTerm( T[kRef] );

      /* Repeat for the right hand side */
      kRef = (nL-1-iSide)*nD + nD-i-1;
      CreepLaw C1right( A[kRef], E[kRef], n[kRef] );
      double A1right = C1right.temperatureTerm( T[kRef] );

      /* Loop through down-dip grid coordinates */
      for( int j=0; j<nD; j++ )
        {
          /* Set the alternative A to mimic what is there */
          CreepLaw C2left( A[i*nD+j], E[i*nD+j], n[kRef] );
          double A2left = C2left.temperatureTerm( T[i*nD+j] );

          CreepLaw C2right( A[iRev*nD+j], E[iRev*nD+j], n[kRef] );
          double A2right = C2right.temperatureTerm( T[iRev*nD+j] );

          /* Choose whichever is weaker */
          A[i*nD+j] = std::max( A1left, A2left )/(A2left/A_amp);
          A[iRev*nD+j] = std::max( A1right, A2right )/(A2right/A_amp);

        } // END Loop through depth

    }

    return;
}

/* similar boundary condition in Ben-Zion 1996 */
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
                   const double &Tgradient)
{
  if( A.size() != nL*nD || E.size() != nL*nD ){
    std::cerr << "side creep ERROR: size of A and E must be nL*nD" << std::endl;
    std::terminate();
  }
  BackgroundTemperatureProfile T( Tsurface, Tgradient);

  double xBDRev = nL*cellLength-xBD;

  // int iSide = floor(xBD/cellLength-0.5);
  int iSide = nL/2;
  for( int i=0; i<iSide; i++ )
    {
      int iRev = nL-i-1; // nL=128, half=64
      
      for( int j=0; j<nD; j++ ){
        double depth = ((double)j + 0.5)*cellHeight;
        double length = ((double)i + 0.5)*cellLength;
        double lengthRev = ((double)iRev + 0.5)*cellLength;
        int iCell_left = i*nD+j;
        int iCell_right = iRev*nD+j;
        // std::cerr<<"check value i,j : "<<i<<", "<<j<<std::endl;
        // std::cerr <<"check value: loading rate: "<<loadingStrainRate<<std::endl;
        // std::cerr <<"check value stress term: "<<pow(sigma_eff,3)<<std::endl;
        // std::cerr <<"check value temp term: "<<exp(-E[iCell_left]/(PhysicalConstants::R_g * T(depth) ))<<std::endl;
        // std::cerr <<"check value T(depth): "<<T(depth)<<std::endl;
        // std::cerr <<"check value E["<<iCell_left<<"]: "<<E[iCell_left]<<std::endl;
        // std::cerr <<"check value constant B: "<<3*log(tau_ratio)/Bterm_length<<std::endl;
        // std::cerr <<"check value std::max: "<<std::max(std::max((depth-zBD),0.0),std::max((length-xBD),0.0))<<std::endl;

        // could change activation energy at shallow depth

        A[iCell_left] = loadingStrainRate/(pow(sigma_eff,3)*exp(-E[iCell_left]/(PhysicalConstants::R_g * T(depth) )))
              *exp(3*log(tau_ratio)/Bterm_length*std::max(std::max((depth-zBD),0.0),std::max((xBD-length),0.0)));

        A[iCell_right] = loadingStrainRate/(pow(sigma_eff,3)*exp(-E[iCell_right]/(PhysicalConstants::R_g * T(depth) )))
              *exp(3*log(tau_ratio)/Bterm_length*std::max(std::max((depth-zBD),0.0),std::max((lengthRev-xBDRev),0.0)));

        // if( j == 0 ){
        //   E[iCell_left] = 0.87*E[iCell_left];
        //   E[iCell_right] = 0.87*E[iCell_right];
        // }
        // if( j == 1){
        //   E[iCell_left] = 0.9*E[iCell_left];
        //   E[iCell_right] = 0.9*E[iCell_right];
        // }
        // if( j == 2){
        //   E[iCell_left] = 0.933*E[iCell_left];
        //   E[iCell_right] = 0.933*E[iCell_right];
        // }
        // if( j ==3){
        //   E[iCell_left] = 0.967*E[iCell_left];
        //   E[iCell_right] = 0.967*E[iCell_right];
        // }
      

      }
    }
}

// side_creep.cpp ends here
