#ifndef STRENGTH_PROFILE_H_
#define STRENGTH_PROFILE_H_

#include <vector>

//------------------------------------------------------------------------------
/**
    Class for defining the static (frictional) strength vs depth
*/
class StaticStrengthProfile{

public:
  /* Constructor *//*
  StaticStrengthProfile(
			const double &cohesion, // frictional cohesion
			const double &coeffFriction, // static frictional coefficient
			const double &dsigmaNdz, // normal stress gradient with depth
			const double &dsigmaPdz ):  // pore presure gradient
    _C(cohesion), _mu(coeffFriction), _dsigmadz( dsigmaNdz - dsigmaPdz){}
  */

  /* Constructor for effective normal stress */
  StaticStrengthProfile(
			const double &cohesion, // frictional cohesion
			const std::vector<double> &coeffFriction, // static frictional coefficient
			const double &dsigmadz, // effective normal stress gradient with depth
      const int &nL, const int &nD,
      const double &cellHeight
		       ):
      _C(cohesion), _mu(coeffFriction), _dsigmadz( dsigmadz ), _nL(nL), _nD(nD), _cellHeight(cellHeight) {}

public:
  void getStaticStrength( std::vector<double> &taus) const;

  double getCellStaticStrength(const int &i, const int &j) const;



 private:
  /* Properties needed to compute the strength profile */
  double _C; // Cohesion
  std::vector<double> _mu;  // static coefficient of friction
  double _dsigmadz; // Effective normal stress gradient with depth
  int _nL;
  int _nD;
  double _cellHeight;

};


#endif /*STRENGTH_PROFILE_H_*/
