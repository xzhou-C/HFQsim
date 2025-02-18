#include "strength_profile.h"

//------------------------------------------------------------------------------
/* Define profile on a fault grid, column ordered */
void StaticStrengthProfile::getStaticStrength(std::vector<double> &taus)
  const
{
  if( taus.size() != _nL * _nD){
    taus.resize( _nL*_nD );
  }
  for(unsigned int i=0; i<_nL; i++){
    for(unsigned int j=0; j<_nD; j++){
      taus[i*_nD+j] = getCellStaticStrength(i,j);
    }
  }
  return;
}

double StaticStrengthProfile::getCellStaticStrength( const int &i, 
                                                     const int &j)
  const
  {
    return _C + _mu[i*_nD+j] * _dsigmadz * _cellHeight * ((double)j + 0.5);
  }

//------------------------------------------------------------------------------

