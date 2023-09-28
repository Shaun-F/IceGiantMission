#ifndef LISA_H_INCLUDED
#define LISA_H_INCLUDED

#include "Constants.hpp"
#include "utils.hpp"
#include "VariablesParameters.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <filesystem> //for capturing absolute directory of compiled executable
using namespace std;







class LISA {
private:
  using TransferData=std::vector<std::vector<high_prec_t>>;

  high_prec_t m_fm{3.168753575e-8}; //LISA modulation frequency
  high_prec_t m_Tobs{}; //observation time in years
  high_prec_t m_Larm{}; //lisa arm length in meters
  int m_NC {}; //number of data channels
  high_prec_t m_ecc { static_cast<high_prec_t>(m_Larm)/(2.0*sqrt(3.0)*Constants::AU) }; //maintain quasi-equilateral triangle configuration
  high_prec_t m_fstar { Constants::C/(2.0*Constants::PI*m_Larm) };
  high_prec_t m_lisaOrbitFrequency { Constants::omegaE }; //omega_E in python script

  TransferData m_transferdata{};
  bool m_NUMERICALTRANSFERDATA_SUCCESS{};

public:
  LISA(high_prec_t Tobs=4.0*Constants::yr, high_prec_t Larm=2.5e9, int NC=2, std::string transferfilename = "source/R.txt");

  void LoadTransfer(std::string file_name, std::vector< std::vector<high_prec_t>>& transferdata); //load numerical transfer function
  high_prec_t EvaluateTransfer(high_prec_t frequency); //evaluate the transfer function at a single point

  high_prec_t strainPSD(high_prec_t freq); //string power spectral density
  high_prec_t snC(high_prec_t freq); //galactic confusion noise
  high_prec_t sN(high_prec_t freq); //sensitivity curve
  high_prec_t pn_WC(high_prec_t freq); //power spectral density with confusion

  std::vector<high_prec_t> patternFunctions(params&);
  std::vector<high_prec_t> patternFunctions(high_prec_t, high_prec_t, high_prec_t); //overloaded
  std::vector<high_prec_t> timedepPatternFunctions(params&, high_prec_t);

};







#endif //LISA_H_INCLUDED
