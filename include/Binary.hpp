#ifndef BINARY_H_INCLUDED
#define BINARY_H_INCLUDED
#include "LISA.hpp"
#include "Constants.hpp"
#include "utils.hpp"
#include "VariablesParameters.hpp"
using namespace std;





class Binary: public LISA{
private:
  high_prec_t m_chirpMass{};
  high_prec_t m_f1{}; //time derivative of phase evolution
  high_prec_t m_a0{};
  params& m_paramStruct;
  std::vector<high_prec_t> m_gwAmplitude{};

public:
  Binary()=delete; //require insertion of params
  //constructor
  Binary(params& params);

  high_prec_t getChirpMass(){ return m_chirpMass;};
  std::vector<high_prec_t> getGWAmplitude(){ return m_gwAmplitude;};

  //methods are static to easily allow the use of GSL integrators.
  std::vector<high_prec_t>  GWAmplitudes(params&);
  high_prec_t FrequencyK(params&);
  high_prec_t Separation(params&);
  high_prec_t Amplitude(params&, high_prec_t);
  high_prec_t polarizationPhase(params&, high_prec_t);
  high_prec_t dopplerPhase(params&, high_prec_t);
  high_prec_t gwFrequency(params&, high_prec_t);
  high_prec_t phase(params&, high_prec_t);
  high_prec_t strain(params&, high_prec_t);
  high_prec_t freqDeriv(params&, ParameterVariables, high_prec_t);
  high_prec_t freqDerivInt(params&, ParameterVariables, high_prec_t);
  high_prec_t dStrain(params&, ParameterVariables, high_prec_t);


  //doppler mission specific functions
  high_prec_t lineOfSightAngle(params&);
  high_prec_t detectorStrain(params&, high_prec_t);
  high_prec_t psiBar(params&, high_prec_t);
  high_prec_t dPsiBar(params&, high_prec_t);

  #ifdef COMPILE_AS_PROGRAM
  
  //generate fisher matrix
  high_prec_t FisherElement(params&, ParameterVariables, ParameterVariables, int);
  //static high_prec_t Integrand(high_prec_t, void*);

  void genFisherMatrix(params&, FisherMatrix&, high_prec_t&, int);

  #endif //COMPILE_AS_PROGRAM


};



#endif //BINARY_H_INCLUDED
