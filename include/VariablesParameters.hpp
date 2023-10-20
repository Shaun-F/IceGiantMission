#ifndef VARIABLESPARAMETERS_H_INCLUDED
#define VARIABLESPARAMETERS_H_INCLUDED
#include <vector>
#include <string>
//Boost libraries for multiprecision
#include <boost/math/constants/constants.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
//using boost::multiprecision::cpp_dec_float_50;
using high_prec_t=double; //boost::multiprecision::cpp_dec_float<25>;
using namespace std;



enum class ParameterVariables{
  K,
  P,
  phiP,
  thetaS,
  phiS,
  thetaL,
  phiL,
  lnA,
  f1,
  f0,
  NUM_Variables,
};

struct analysisParams{//struct containing analysis parameters
  std::vector<high_prec_t> PeriodList {0.01,0.03,0.05,0.07,0.09,0.1,0.3,0.5,0.7,0.9,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,20.,30.,40.,50.,60.,70.,80.,90. }; //list of CBP periods to iterate over
};



struct params{ //struct containing parameters for specific system
  high_prec_t thetaL{0.};
  high_prec_t phiL{0.};
  high_prec_t thetaS{0.};
  high_prec_t phiS{0.};
  high_prec_t M1{0.}; //kg
  high_prec_t M2{0.}; //kg
  high_prec_t MP{0.}; //kg
  high_prec_t P{0.}; //seconds
  high_prec_t Tobs{4.0*Constants::yr}; //seconds
  high_prec_t Larm{2.5e9}; //meters
  high_prec_t thetaP{0.};
  high_prec_t phiP{0.};
  high_prec_t freqGW{0.}; //hz
  high_prec_t sourceDistance{0.}; //meters
  high_prec_t cosIota{0.0};
  high_prec_t gravAmplitude{0.0};
  high_prec_t DerivativeDelta{1e-6};
  high_prec_t chirpMass{};
  high_prec_t f1{};
  high_prec_t K{};
  high_prec_t ig_direction{};
  high_prec_t LineOfSightAngle{};
  high_prec_t lightTwoWayTime{}; //two way light travel time for doppler mission
  double PSDlevel{0};
  double allanDeviation{0};
  double relativeSNR{1e-1};
  int periodSamples{35};
  int NC{2};
  int LISAAlpha{0};  //denotes which beam pattern function were looking at (1, 2, or 3(both))
  char mode{'s'};
  std::vector<high_prec_t> gwAmplitude{};
  std::string mission{"LISA"}; //options: ["LISA", "IceGiant"] do we look at the LISA mission or Ice Giant mission

};






#endif //VARIABLESPARAMETERS_H_INCLUDED
