#include "Constants.hpp"
#include "LISA.hpp"
#include "VariablesParameters.hpp"
#include "VariablesParameters.hpp"
#include <string>
#include <vector>
#include <fstream> //for reading files
#include <iostream> //for writing to console
#include <filesystem>
using namespace std;


//constructor for LISA class
LISA::LISA(high_prec_t Tobs, high_prec_t Larm, int NC, std::string transferfilename): m_Tobs{Tobs}, m_Larm{Larm}, m_NC{NC} {
    //calculate the transfer data and store in private member variable
    auto dataPath {std::filesystem::current_path()}; //get current working directory
    dataPath /= transferfilename;
    LoadTransfer(dataPath, m_transferdata);
};

void LISA::LoadTransfer(std::string file_name, std::vector<std::vector<high_prec_t>>& transferFunctionData){

  std::ifstream TransferFile{file_name};

  if (!TransferFile) { //if transfer file fails to open, set flag
    std::cout << "Warning: Transfer data failed to load at " << file_name << " . Setting NUMERICALTRANSFERDATA_SUCCESS flag to false and using analytic approximations." <<std::endl;
    m_NUMERICALTRANSFERDATA_SUCCESS = false; //failed to load file
  } else { //file successfully opened
    std::cout << "Successfully loaded Transfer Data"<<std::endl;
    m_NUMERICALTRANSFERDATA_SUCCESS = true; //successfully loaded file


    high_prec_t freqPoint{}; //independent variable
    high_prec_t valuePoint{}; //dependent variable
    while (!TransferFile.eof()) { //while we're still scanning the file
      TransferFile >> freqPoint >> valuePoint; // parse line to the two variables
      std::vector<high_prec_t> tmpPoint {freqPoint*m_fstar, valuePoint*m_NC }; //multiply by fstar to get frequency, multiply by NC to get improvements due to more data channels
      transferFunctionData.push_back(tmpPoint);//add data to container
    }; //while
  }; //ifelse
};


high_prec_t LISA::EvaluateTransfer(high_prec_t frequency){

  high_prec_t returnValue{};
  if (m_NUMERICALTRANSFERDATA_SUCCESS){ //if data file is provided
    for (int iter{0}; iter<m_transferdata.size(); ++iter){
      if (m_transferdata[iter][0]>frequency){

        //linear interpolation
        high_prec_t lowerFrequency { m_transferdata[iter-1][0]};
        high_prec_t higherFrequency { m_transferdata[iter][0]};
        high_prec_t lowerValue { m_transferdata[iter-1][1]};
        high_prec_t higherValue { m_transferdata[iter][1]};


        high_prec_t slope { (higherValue - lowerValue)/(higherFrequency - lowerFrequency) };
        high_prec_t intercept { (lowerValue*higherFrequency - higherValue*lowerFrequency)/(higherFrequency - lowerFrequency) };

        returnValue = slope*frequency + intercept;
        break;
      }; //inner if
    }; //for loop
  } else {
    returnValue = (3./20.)*(1./(1. + (6./10.)*(frequency*m_fstar)*(frequency*m_fstar)))*m_NC;
  };

  return returnValue;
};


high_prec_t LISA::strainPSD(high_prec_t freq){
  high_prec_t P_oms { pow(1.5e-11,2.) * (1. + pow(2.0e-3/freq,4)) };

  high_prec_t P_acc_t1 {pow(3.0e-15,2)};
  high_prec_t P_acc_t2 {1. + pow((0.4e-3/freq),2)};
  high_prec_t P_acc_t3 { 1. + pow(freq/(8.0e-3),4)};
  high_prec_t P_acc {P_acc_t1*P_acc_t2*P_acc_t3};

  high_prec_t Pn{ (P_oms + 2.*(1. + pow(cos(freq/m_fstar),2)) * P_acc/(pow(2.*Constants::PI*freq,4)))/(pow(m_Larm,2)) };

  return Pn;
}


high_prec_t LISA::snC(high_prec_t freq){
  high_prec_t alpha{};
  high_prec_t beta{};
  high_prec_t gamma{};
  high_prec_t kappa{};
  high_prec_t f_knee{};
  if (m_Tobs<0.75*Constants::yr){
    alpha  = 0.133;
    beta   = 243.;
    kappa  = 482.;
    gamma  = 917.;
    f_knee = 2.58e-3;
  } else if (m_Tobs>0.75*Constants::yr && m_Tobs<1.5*Constants::yr) {
    alpha  = 0.171;
    beta   = 292.;
    kappa  = 1020.;
    gamma  = 1680.;
    f_knee = 2.15e-3;
  } else if (m_Tobs>1.5*Constants::yr && m_Tobs<3.0*Constants::yr) {
    alpha  = 0.165;
    beta   = 299.;
    kappa  = 611.;
    gamma  = 1340.;
    f_knee = 1.73e-3;
  } else {
    alpha  = 0.138;
    beta   = -221.;
    kappa  = 521.;
    gamma  = 1680.;
    f_knee = 1.13e-3;
  };

  high_prec_t Amp {1.8e-44 / m_NC};
  high_prec_t Sc { 1. + tanh(gamma*(f_knee - freq))};
  Sc *= exp(-pow(freq, alpha) + beta*freq*sin(kappa*freq));
  Sc *= Amp*pow(freq, -7./3.);

  return Sc;
}

high_prec_t LISA::sN(high_prec_t freq){
  //LISA strain sensitivity curve
  high_prec_t TransferValue{EvaluateTransfer(freq)};

  high_prec_t Sn {strainPSD(freq)/TransferValue + snC(freq)};
  return Sn;
}


high_prec_t LISA::pn_WC(high_prec_t freq){


  high_prec_t TransferValue{EvaluateTransfer(freq)};
  high_prec_t PnC { strainPSD(freq) + snC(freq)*TransferValue };

  return PnC;
}


std::vector<high_prec_t> LISA::patternFunctions(params& params){
  high_prec_t thetaSv = params.thetaS;
  high_prec_t phiSv = params.phiS;
  high_prec_t thetaL { params.thetaL };
  high_prec_t phiL { params.phiL };
  high_prec_t cosIota { params.cosIota };
  high_prec_t psiSv = psiS(thetaL, thetaSv, phiL, phiSv, cosIota, 0);

  std::vector<high_prec_t> patternFuncs {};
  high_prec_t termA { 0.5 * (1. + pow(cos(thetaSv),2)) * cos(2*phiSv)};
  high_prec_t termB { cos(thetaSv) * sin(2*phiSv) };

  patternFuncs = {termA * cos(2*psiSv) - termB * sin(2*psiSv),
                  termA * sin(2*psiSv) + termB * cos(2*psiSv)};
  return patternFuncs;
}

//overload LISA::patternFunctions method
std::vector<high_prec_t> LISA::patternFunctions(high_prec_t thetaS, high_prec_t phiS, high_prec_t psiS){

  std::vector<high_prec_t> patternFuncs {};
  high_prec_t termA { 0.5 * (1. + pow(cos(thetaS),2)) * cos(2*phiS)};
  high_prec_t termB { cos(thetaS) * sin(2*phiS) };

  patternFuncs = {termA * cos(2*psiS) - termB * sin(2*psiS),
                  termA * sin(2*psiS) + termB * cos(2*psiS)};
  return patternFuncs;
}


std::vector<high_prec_t> LISA::timedepPatternFunctions(params& params, high_prec_t time){
  if (params.mission=="LISA"){
    high_prec_t thetaL = params.thetaL;
    high_prec_t phiL = params.phiL;
    high_prec_t thetaSv = params.thetaS;
    high_prec_t phiSv = params.phiS;
    high_prec_t cosIota = CosIota(params);
    int LISAAlpha = params.LISAAlpha;

    high_prec_t thetaSval { acos(cosThetaS(thetaSv, phiSv, time)) };
    high_prec_t phiSval { phiS(thetaSv, phiSv, time) - Constants::PI/4 * (LISAAlpha-1) };
    high_prec_t psiSval { psiS(params, time) };
    return patternFunctions(thetaSval, phiSval, psiSval);

  } else if (params.mission == "IceGiant") {
    high_prec_t thetaSv = params.thetaS;
    high_prec_t phiSv = params.phiS;
    high_prec_t psiSv = psiS(params, time);
    high_prec_t igDirect = params.ig_direction;
    std::vector<high_prec_t> BPF { cos(thetaSv)*cos(thetaSv)*cos(phiSv - igDirect)*cos(phiSv - igDirect) - sin(phiSv - igDirect)*sin(phiSv - igDirect), 2*cos(thetaSv)*sin(phiSv - igDirect)*cos(phiSv-igDirect)};
    std::vector<high_prec_t> result { BPF[0]*cos(2*psiSv) - BPF[1]*sin(2*psiSv),
                                 BPF[0]*sin(2*psiSv) + BPF[1]*cos(2*psiSv)   };
    return result;
  } else {
    assert(("Mission parameter not known. ", false));
    return {0.,0.};
  }

}
