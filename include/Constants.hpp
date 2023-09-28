#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED

using Constant=double;
namespace Constants{
  inline constexpr Constant PI {3.141592653589793}; //PI
  inline constexpr Constant AU {149597870700.0}; //AU's to meters
  inline constexpr Constant pc {3.085677581491367e+16};//parsec to meters
  inline constexpr Constant yr {60*60*24*365.25}; //seconds in year
  inline constexpr Constant day {60*60*24}; //second in day
  inline constexpr Constant MJ {1.8981245973360505e+27}; //mass of jupiter in kg
  inline constexpr Constant MS {1.9884098e30}; //mass of sun in kg
  inline constexpr Constant Gnewt {6.6743e-11}; //gravitational constant
  inline constexpr Constant C {299792458.0}; //speed of light
  inline constexpr Constant rS {2*Gnewt/(C*C)}; //schwarzschild radius
  inline constexpr Constant omegaE {2*PI/yr}; //earth orbit frequency
};
#endif //CONSTANTS_H_INCLUDED
