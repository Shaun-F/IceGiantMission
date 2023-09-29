#include "Binary.hpp"
#include "Constants.hpp"
#include "LISA.hpp"
#include "VariablesParameters.hpp"
#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <string>
#include <omp.h>
#include <cassert>
using namespace std;


Binary::Binary(params& params): m_paramStruct(params), LISA(params.Tobs, params.Larm, params.NC){
  m_chirpMass = pow( params.M1 * params.M2, 3./5. ) / pow( params.M1 + params.M2 , 1./5. );
  m_f1 = 96./5. * pow(Constants::PI, 8./3.) * pow( params.freqGW , 11./3.) * pow(Constants::Gnewt * m_chirpMass/pow(Constants::C,3.), 5./3.);

  params.chirpMass = m_chirpMass;
  params.f1 = m_f1;
  params.a0 = m_a0;
}

std::vector<high_prec_t> Binary::GWAmplitudes(params& params){
  //!!!!!!!!!!!!!!!! second element should be negative
  high_prec_t preFactor { (4.0/params.sourceDistance) * pow(Constants::Gnewt*params.chirpMass / (pow(Constants::C,2)), 5./3.) * pow( Constants::PI*params.freqGW/Constants::C, 2./3.) };
  return {preFactor * (1. + CosIota(params)*CosIota(params))/2., -preFactor*CosIota(params)};
}

high_prec_t Binary::FrequencyK(params& params){
  return params.MP*sin(params.thetaP)/Constants::C * pow( 2*Constants::PI*Constants::Gnewt/params.P , 1./3.) * pow(1./(params.M1 + params.M2 + params.MP), 2./3.);
}



high_prec_t Binary::Separation(params& params){
  high_prec_t value { Constants::Gnewt * (params.M1 + params.M2) / (pow(Constants::PI * params.freqGW,2)) };
  return pow(value,1/6);
}

high_prec_t Binary::Amplitude(params& params, high_prec_t time){
  std::vector<high_prec_t> lisaTimeDepPatterFuncs { timedepPatternFunctions(params,time)};
  high_prec_t norm { sqrt( GWAmplitudes(params)[0]*GWAmplitudes(params)[0]*lisaTimeDepPatterFuncs[0]*lisaTimeDepPatterFuncs[0] + GWAmplitudes(params)[1]*GWAmplitudes(params)[1]*lisaTimeDepPatterFuncs[1]*lisaTimeDepPatterFuncs[1])};
  return norm;
}

high_prec_t Binary::polarizationPhase(params& params, high_prec_t time){
  std::vector<high_prec_t> lisaTimeDepPatterFuncs { timedepPatternFunctions(params,time)};
  std::vector<high_prec_t> tmpVec { GWAmplitudes(params)[0]*lisaTimeDepPatterFuncs[0], GWAmplitudes(params)[1]*lisaTimeDepPatterFuncs[1] };
  return atan2(-tmpVec[1], tmpVec[0]);
}

high_prec_t Binary::gwFrequency(params& params, high_prec_t time){
  //frequency of gravitational wave including perturbation due to circumbinary planet
  return (params.freqGW + params.f1*time) * (1 - FrequencyK(params) * cos(2.*Constants::PI*time/params.P + params.phiP));
}

high_prec_t Binary::dopplerPhase(params& params, high_prec_t time){
  high_prec_t freq = gwFrequency(params,time);
  return 2*Constants::PI * freq / Constants::C * Constants::AU * sin(params.thetaS) * cos(Constants::omegaE * time - params.phiS);
}

high_prec_t Binary::phase(params& params, high_prec_t time){
  //analytic integral from 0 to t of GW frequency
  high_prec_t freq = gwFrequency(params,time);
  high_prec_t T1 { 2*time*(2*freq + params.f1*time) };
  high_prec_t T2 {  params.f1 * params.P * cos(params.phiP) };
  high_prec_t T3 { -params.f1 * params.P * cos(2.*Constants::PI*time/params.P + params.phiP) };
  high_prec_t T4 { 2 * freq * Constants::PI * sin(params.phiP) };
  high_prec_t T5 { -2*Constants::PI * (freq + params.f1*time) *sin(2.*Constants::PI*time/params.P + params.phiP) };
  return 1./4. * (T1 + (FrequencyK(params) * params.P)/pow(Constants::PI, 2) * (T2 + T3 + T4 + T5));

}

high_prec_t Binary::lineOfSightAngle(params& params){
  return sin( params.thetaS ) * ( cos(params.ig_direction) * cos(params.phiS) + sin(params.ig_direction)*sin(params.phiS));

}


high_prec_t Binary::detectorStrain(params& params, high_prec_t time){
  assert(params.mission=="IceGiant");
  return Amplitude(params,time)*cos(2*Constants::PI*phase(params,time) + polarizationPhase(params,time));
}

high_prec_t Binary::strain(params& params, high_prec_t time){
  high_prec_t val {sqrt(3)/2.};
  if (params.mission=="LISA"){
    return val*Amplitude(params, time)*cos(2.*Constants::PI*phase(params, time) + polarizationPhase(params, time) + dopplerPhase(params, time));
  } else if (params.mission == "IceGiant") {
    return (lineOfSightAngle(params)-1)/2 * psiBar(params,time) - lineOfSightAngle(params)*psiBar(params,time-(1+lineOfSightAngle(params))/2 * params.lightTwoWayTime) + (1+lineOfSightAngle(params))/2*psiBar(params,time-params.lightTwoWayTime);
  };
  return 0.;
}

high_prec_t Binary::psiBar(params& params, high_prec_t time){
  assert(params.mission=="IceGiant"); //ERROR: psiBar not defined for missions other than IceGiant
  return detectorStrain(params,time)/(1 - lineOfSightAngle(params)*lineOfSightAngle(params));
}

high_prec_t Binary::dPsiBar(params& params, high_prec_t time){
  assert(params.mission=="IceGiant"); //ERROR: dPsiBar not defined for missions other than IceGiant
  return -Amplitude(params,time)*sin(2*Constants::PI*phase(params,time) + polarizationPhase(params,time))/(1-lineOfSightAngle(params)*lineOfSightAngle(params));

}


high_prec_t Binary::freqDeriv(params& params, ParameterVariables myVar, high_prec_t time){
  //calculate derivative of frequency with respect to parameters specified by the enum ParameterVariables
  high_prec_t result{};
  switch(myVar){
    case ParameterVariables::K:
      result = -(params.freqGW + params.f1*time)*cos(2.*Constants::PI*time/params.P + params.phiP);
      break;
    case ParameterVariables::P:
      result = -1./(params.P*params.P) * 2*FrequencyK(params)*Constants::PI*time*(params.freqGW + params.f1*time)*sin(2.*Constants::PI*time/params.P + params.phiP);
      break;
    case ParameterVariables::phiP:
      result = FrequencyK(params)*(params.freqGW + params.f1*time)*sin(2.*Constants::PI*time/params.P + params.phiP);
      break;
    case ParameterVariables::f0:
      result = 1 - FrequencyK(params)*cos(2.*Constants::PI*time/params.P + params.phiP);
      break;
    case ParameterVariables::f1:
      result = time*(1 - FrequencyK(params)*cos(2.*Constants::PI*time/params.P + params.phiP));
      break;
    default:
      std::cout << "ERROR: variable not valid" << std::endl;
      result = 0.0;
      break;
    };
  return result;
}

high_prec_t Binary::freqDerivInt(params& params, ParameterVariables myVar, high_prec_t time){
  //calculate derivative of frequency with respect to parameters specified by the enum ParameterVariables
  high_prec_t result{};
  switch(myVar){
    case ParameterVariables::K:
      result = (params.P/(4.*Constants::PI*Constants::PI)) * (params.f1 * params.P * cos(params.phiP) - params.f1 * params.P*cos(2.*Constants::PI*time/params.P + params.phiP) + 2*params.freqGW*Constants::PI*sin(params.phiP) - 2*Constants::PI*(params.freqGW + params.f1*time) * sin(2.*Constants::PI*time/params.P + params.phiP));
      break;
    case ParameterVariables::P:
      result = (FrequencyK(params)/(2.*params.P*Constants::PI*Constants::PI)) * (params.f1*params.P*params.P *cos(params.phiP) + (-params.f1*params.P*params.P + 2*params.freqGW*Constants::PI*Constants::PI*time + 2*params.f1*Constants::PI*Constants::PI*time*time)*cos(2.*Constants::PI*time/params.P + params.phiP) + params.P*Constants::PI*(params.freqGW*sin(params.phiP) - (params.freqGW + 2*params.f1*time)*sin(2.*Constants::PI*time/params.P + params.phiP)) );
      break;
    case ParameterVariables::phiP:
      result = (FrequencyK(params)*params.P/(4.*Constants::PI*Constants::PI)) * (2*params.freqGW*Constants::PI*cos(params.phiP) - 2*Constants::PI*(params.freqGW + params.f1*time)*cos(2.*Constants::PI*time/params.P + params.phiP) + params.f1*params.P*(-sin(params.phiP) + sin(2.*Constants::PI*time/params.P + params.phiP)));
      break;
    case ParameterVariables::f0:
      result = time - (1./Constants::PI)*params.P*FrequencyK(params)*cos(2.*Constants::PI*time/params.P + params.phiP)*sin(Constants::PI*time/params.P);
      break;
    case ParameterVariables::f1:
      result = (1./2.)*time*time - (FrequencyK(params)*params.P/(4.*Constants::PI*Constants::PI))*(-params.P*cos(params.phiP) + params.P*cos(2.*Constants::PI*time/params.P + params.phiP) + 2*Constants::PI*time*sin(2.*Constants::PI*time/params.P + params.phiP));
      break;
    default:
      std::cout << "ERROR: variable not valid" << std::endl;
      result = 0.0;
      break;
    };
  return result;
}


high_prec_t Binary::dStrain(params& myParams, ParameterVariables var, high_prec_t time){
  high_prec_t result {};
  high_prec_t val {sqrt(3.)/2.};

  switch (var){
    case ParameterVariables::K:
    case ParameterVariables::P:
    case ParameterVariables::phiP:
    case ParameterVariables::f1:
    case ParameterVariables::f0:
      if (myParams.mission=="LISA"){
        result = -val*Amplitude(myParams, time)*(2*Constants::PI*freqDerivInt(myParams,var,time) + 2*Constants::PI*freqDeriv(myParams,var,time)/Constants::C * Constants::AU * sin(myParams.thetaS) * cos(Constants::omegaE*time - myParams.phiS))*sin(2*Constants::PI*phase(myParams,time) + polarizationPhase(myParams, time) + dopplerPhase(myParams,time));
      } else if (myParams.mission=="IceGiant") {
        result = (lineOfSightAngle(myParams) - 1)/2*2*Constants::PI*freqDerivInt(myParams,var,time)*dPsiBar(myParams,time) - lineOfSightAngle(myParams)*2*Constants::PI*freqDerivInt(myParams,var,time - (lineOfSightAngle(myParams)+1)/2*myParams.lightTwoWayTime)*dPsiBar(myParams,time - (lineOfSightAngle(myParams)+1)/2*myParams.lightTwoWayTime) + (lineOfSightAngle(myParams) + 1)/2*2*Constants::PI*freqDerivInt(myParams,var, time-myParams.lightTwoWayTime)*dPsiBar(myParams,time-myParams.lightTwoWayTime);
      };
      break;
    case ParameterVariables::thetaS:
      { //create scope to hold temporary param copy
      params myTempParams { myParams }; //create local copy of params
      int DerivativeDirection {1};
      if(myParams.thetaS+myParams.DerivativeDelta>Constants::PI){
        DerivativeDirection*=-1; //if ThetaS near Pi, change to left sided derivative to prevent undefined behavior
      };
      myTempParams.thetaS += DerivativeDirection*myParams.DerivativeDelta; //left sided derivative
      result = (strain(myTempParams,time) - strain(myParams, time))/(DerivativeDirection*myParams.DerivativeDelta);
      } //end of temporary scope
      break;
    case ParameterVariables::thetaL:
      { //create scope to hold temporary param copy
      params myTempParams = myParams; //create local copy of params
      int DerivativeDirection {1};
      if(myParams.thetaL+myParams.DerivativeDelta>Constants::PI){
        DerivativeDirection*=-1; //if ThetaS near Pi, change to left sided derivative to prevent undefined behavior
      };
      myTempParams.thetaL += DerivativeDirection*myParams.DerivativeDelta; //left sided derivative
      result = (strain(myTempParams,time) - strain(myParams, time))/(DerivativeDirection*myParams.DerivativeDelta);
      } //end of temporary scope
      break;
    case ParameterVariables::phiS:
      { //create scope to hold temporary param copy
      params myTempParams = myParams; //create local copy of params
      int DerivativeDirection {1};
      if(myParams.phiS+myParams.DerivativeDelta>Constants::PI){
        DerivativeDirection*=-1; //if ThetaS near Pi, change to left sided derivative to prevent undefined behavior
      };
      myTempParams.phiS += DerivativeDirection*myParams.DerivativeDelta; //left sided derivative
      result = (strain(myTempParams,time) - strain(myParams, time))/(DerivativeDirection*myParams.DerivativeDelta);
      } //end of temporary scope
      break;
    case ParameterVariables::phiL:
      { //create scope to hold temporary param copy
      params myTempParams = myParams; //create local copy of params
      int DerivativeDirection {1};
      if(myParams.phiL+myParams.DerivativeDelta>Constants::PI){
        DerivativeDirection*=-1; //if ThetaS near Pi, change to left sided derivative to prevent undefined behavior
      };
      myTempParams.phiL += DerivativeDirection*myParams.DerivativeDelta; //left sided derivative
      result = (strain(myTempParams,time) - strain(myParams, time))/(DerivativeDirection*myParams.DerivativeDelta);
      } //end of temporary scope
      break;
    case ParameterVariables::lnA:
      result = strain(myParams, time);
      break;
    default:
      std::cout << "ERROR in dStrain. Specified variable not recognized." <<std::endl;
      break;
  };

  return result;
}

high_prec_t Binary::FisherElement(params& myParams, ParameterVariables var1, ParameterVariables var2, int numThreads=1){
  high_prec_t totalResult {0}; //result container, zero initialized
  int periodSamples { myParams.periodSamples };
  auto integrand {
		  [&](high_prec_t time)->high_prec_t {
			  return this->dStrain(myParams, var1, time) * this->dStrain(myParams, var2, time); 
		  }
	  };
  //perform integration
  auto function1{
    [&](high_prec_t time) -> high_prec_t {
      return this->dStrain(myParams, var1, time);
    }
  };
  auto function2{
    [&](high_prec_t time) -> high_prec_t {
      return this->dStrain(myParams, var2, time);
    }
  };

  if (myParams.mission=="LISA"){
    high_prec_t Nsamples { myParams.Tobs * myParams.freqGW * periodSamples };
    high_prec_t dx { (high_prec_t)myParams.Tobs/Nsamples};


    /////////////////first beam pattern function
    myParams.LISAAlpha = 1; 
    { //create a code block to localize variables. temperary variables like input1/2, ffts get destroyed at end of codeblock
      //compute input arrays
      std::vector<high_prec_t> input1 { sampleFunction(function1, 0, myParams.Tobs, Nsamples, numThreads) };
      std::vector<high_prec_t> input2 { sampleFunction(function2, 0, myParams.Tobs, Nsamples, numThreads) };

      //caluclate ffts of input arrays
      std::vector<std::vector<high_prec_t>> fft1 { fourierTransform(input1, numThreads) };
      std::vector<std::vector<high_prec_t>> fft2 { fourierTransform(input2, numThreads) };

      //compute integrand by multiplying ffts together
      std::vector<high_prec_t> realPart{};
      std::vector<high_prec_t> complexPart{};
      for (int i{0}; i<Nsamples; ++i){
        realPart.push_back(fft1[i][0]*fft2[i][0] + fft1[i][1]*fft2[i][1]);
        complexPart.push_back(fft1[i][1]*fft2[i][0] - fft1[i][0]*fft2[i][1]);
      };

      //compute integral
      high_prec_t integral1 { myTrapSum(realPart, dx, numThreads) };
      high_prec_t integral1c { myTrapSum(complexPart, dx, numThreads) };
      if (integral1c>1e-10){
        std::cout << "Complex part of integral: " << integral1c << std::endl;
      };

      totalResult = 2*(integral1)/sN(myParams.freqGW);
    };


    /////////////////second beam pattern function
    myParams.LISAAlpha = 2; 
    { //create a code block to localize variables. temperary variables like input1/2, ffts get destroyed at end of codeblock
      //compute input arrays
      std::vector<high_prec_t> input21 { sampleFunction(function1, 0, myParams.Tobs, Nsamples, numThreads) };
      std::vector<high_prec_t> input22 { sampleFunction(function2, 0, myParams.Tobs, Nsamples, numThreads) };

      //caluclate ffts
      std::vector<std::vector<high_prec_t>> fft21 { fourierTransform(input21, numThreads) };
      std::vector<std::vector<high_prec_t>> fft22 { fourierTransform(input22, numThreads) };

      //compute integrand
      std::vector<high_prec_t> realPart2{};
      std::vector<high_prec_t> complexPart2{};
      for (int i{0}; i<Nsamples; ++i){
        realPart2.push_back(fft21[i][0]*fft22[i][0] + fft21[i][1]*fft22[i][1]);
        complexPart2.push_back(fft21[i][1]*fft22[i][0] - fft21[i][0]*fft22[i][1]);
      }

      //compute integral
      high_prec_t integral2 { myTrapSum(realPart2, dx, numThreads) };
      high_prec_t integral2c { myTrapSum(complexPart2, dx, numThreads) };
      if (integral2c>1e-10){
        std::cout << "Complex part of integral2: " << integral2c << std::endl;
      }

      totalResult += 2*(integral2)/sN(myParams.freqGW);
    };

  } else if (myParams.mission=="IceGiant"){
    //array of start values in units of years
    std::array<high_prec_t, 10> startTimes { 0., 10./9., 20./9., 30./9., 40./9., 50./9., 60./9., 70./9., 80./9., 90./9. };
    high_prec_t Nsamples { 40.*Constants::day * myParams.freqGW * periodSamples };
    
    //iterate over start times. 10 40-day observation windows
    for (auto& stime: startTimes){
      stime *= Constants::yr; //convert to seconds
      high_prec_t Integral { myTrapIntegral(integrand, stime, stime + 40.*Constants::day, Nsamples, numThreads) };
      totalResult += Integral;
    };

  } else {
    assert(("ERROR: Mission parameter not recognized", false));
  };

  return totalResult;
}


void Binary::genFisherMatrix(params& myParams, FisherMatrix& out, high_prec_t &LISA_SNR, int numOfThreads=0){
  //parallely calculate elements of fisher matrix

  int maxThreads{};
  if (numOfThreads==0){
    maxThreads = (int)(omp_get_max_threads());
  } else {
    maxThreads = numOfThreads;
  };

  int numCalls { (int)(out.Rows*(out.Rows+1)/2) };
  int numThreads { (int)std::fmin(numCalls, maxThreads) };
  int numNestedThreads { (maxThreads - (maxThreads % numThreads))/numThreads };//use integer division to effectively apply floor function
  if (numNestedThreads*numThreads != omp_get_max_threads() ){
    numNestedThreads = (int)omp_get_max_threads()/numThreads;
  }
  std::cout << "Number outer threads: " << numThreads <<std::endl;
  std::cout << "Number nested threads: " << numNestedThreads <<std::endl;

  //Pre-compute LISA Snr
  { //localize tempParams to code block. Gets destroyed after codeblock
    params tempParams = myParams;
    tempParams.mission="LISA";
    high_prec_t Nsamples { tempParams.Tobs * tempParams.freqGW * tempParams.periodSamples };
    LISA_SNR = std::sqrt((2./sN(myParams.freqGW))*myTrapIntegral( [&](high_prec_t t)->high_prec_t {
                                                  return this->strain(tempParams, t) * this->strain(tempParams,t);
                                                },
                                              0, tempParams.Tobs, Nsamples, maxThreads));
  };

  //if executing IceGiant fisher matrix, compute int h(t) h(t) dt
  high_prec_t IceGiantWaveformNorm;
  if (myParams.mission=="IceGiant"){
    std::array<high_prec_t, 10> startTimes { 0., 10./9., 20./9., 30./9., 40./9., 50./9., 60./9., 70./9., 80./9., 90./9. };
    high_prec_t Nsamples { 40.*Constants::day * myParams.freqGW * myParams.periodSamples };
    for (auto& stime: startTimes){
      stime *= Constants::yr; //convert to seconds
      IceGiantWaveformNorm = myTrapIntegral( [&](high_prec_t t)->high_prec_t {
                                                  return this->strain(myParams, t) * this->strain(myParams,t);
                                                },
                                              stime, stime + 40.*Constants::day, Nsamples, maxThreads);
    };
  };

	omp_set_nested(true); //allow nested parallelism. Important for parallelizing the individual integrations too
  #pragma omp parallel for num_threads(numThreads) collapse(2)
  for(unsigned int i=0; i<out.Rows; ++i){
    for (unsigned int j=i; j<out.Cols; ++j){
      std::string printString { "Calculating Fisher matrix element for parameters [" + paramVarToString(static_cast<ParameterVariables>(i)) + ", " + paramVarToString(static_cast<ParameterVariables>(j)) + "].\n" };
      std::cout << printString;

      //generate result
      high_prec_t result{ FisherElement(myParams, static_cast<ParameterVariables>(i), static_cast<ParameterVariables>(j), numNestedThreads) };
      high_prec_t IntegrationResult;
      if (myParams.mission=="LISA"){
       IntegrationResult = result;
      } else if (myParams.mission=="IceGiant"){
        IntegrationResult = (myParams.relativeSNR*LISA_SNR)*(myParams.relativeSNR*LISA_SNR)*result/IceGiantWaveformNorm;
      } else {
        assert(("Mission parameter not recognized", false));
      }
      
      //add result and errors to output matrix and use symmetry of Fisher matrix to automatically populate symmetric indices
      out(i,j) = IntegrationResult;
      if (i!=j){
        out(j,i) = out(i,j);
      };
    };
  };
}


/*

high_prec_t Binary::FisherElement(params& myParams, ParameterVariables var1, ParameterVariables var2, int numThreads=1){
  high_prec_t totalResult; //result container
  int periodSamples { myParams.periodSamples };
  //perform integration
  auto integrand {
		  [&](high_prec_t time)->high_prec_t {
			  return this->dStrain(myParams, var1, time) * this->dStrain(myParams, var2, time); 
		  }
	  };

  if (myParams.mission=="LISA"){
    high_prec_t Nsamples { myParams.Tobs * myParams.freqGW * periodSamples };

    myParams.LISAAlpha = 1; //first beam pattern function
    high_prec_t Integral1 { myTrapIntegral(integrand, 0., myParams.Tobs, Nsamples, numThreads) };
    
    myParams.LISAAlpha = 2; //second beam pattern function
    high_prec_t Integral2 { myTrapIntegral(integrand, 0., myParams.Tobs, Nsamples, numThreads) };

    totalResult = 2*(Integral1 + Integral2)/sN(myParams.freqGW);

  } else if (myParams.mission=="IceGiant"){
    //array of start values in units of years
    std::array<high_prec_t, 10> startTimes { 0., 10./9., 20./9., 30./9., 40./9., 50./9., 60./9., 70./9., 80./9., 90./9. };

    high_prec_t Nsamples { 40.*Constants::day * myParams.freqGW * periodSamples };
    
    //iterate over start times. 10 40-day observation windows
    for (auto& stime: startTimes){
      stime *= Constants::yr; //convert to seconds
      high_prec_t Integral { myTrapIntegral(integrand, stime, stime + 40.*Constants::day, Nsamples, numThreads) };
      totalResult += 2*Integral/1.;
    };
  } else {
    assert(("ERROR: Mission parameter not recognized", false));
  };

  return totalResult;
}
*/

/*
struct IntegrandParams{
  IntegrandParams()=default; //default constructor
  params myParams;
  ParameterVariables var1;
  ParameterVariables var2;
  Binary& BinaryObject;
};
high_prec_t Integrand(high_prec_t time, void* paramsFromAbove){
    IntegrandParams& _inner_params = *(IntegrandParams*)paramsFromAbove;
    high_prec_t result {_inner_params.BinaryObject.dStrain(_inner_params.myParams, _inner_params.var1,time)*_inner_params.BinaryObject.dStrain(_inner_params.myParams, _inner_params.var2,time)};
    return result;
};


std::vector<high_prec_t> Binary::FisherElement(params& myParams, ParameterVariables var1, ParameterVariables var2){

  high_prec_t result1, error1, result2, error2;
  //setup integration memory space
  gsl_integration_workspace* workSpace {gsl_integration_workspace_alloc(1000000)};
  gsl_function F1;
  gsl_function F2;
  IntegrandParams customParams{myParams, var1, var2, *this};
  F1.function = Integrand;
  F1.params = &customParams;

  high_prec_t totalResult{0.};
  high_prec_t totalError{0.};

  high_prec_t epsilonRelative{1e-2};
  high_prec_t epsilonAbsolute{0.};
  int status1{};
  int status2{};

  //QAG adaptive integration

  //turn off default error handling
  gsl_set_error_handler_off();

  //perform integration
  if (myParams.mission=="LISA"){
    customParams.myParams.LISAAlpha = 1; //first beam pattern function
    status1 = gsl_integration_qag(&F1, 0., myParams.Tobs, epsilonAbsolute, epsilonRelative, 1000000, 6, workSpace, &result1, &error1);
    customParams.myParams.LISAAlpha = 2; //second beam pattern function
    status2 = gsl_integration_qag(&F1, 0., myParams.Tobs, epsilonAbsolute, epsilonRelative, 1000000, 6, workSpace, &result2, &error2);

    totalResult = 2*(result1 + result2)/sN(myParams.freqGW);
    totalError = sqrt(error1*error1 + error2*error2);  //error handling

    //error handling
    if (status1==GSL_EROUND || status2==GSL_EROUND){
      std::cout << "Round off error detected. \n result1 = " << result1 << "  result2 = " << result2 << "  error1 = " << error1 << "  error2 = " << error2 << "  relative error1 = " << error1/result1 << "  relative error2 = " << error2/result2 << std::endl;
      std::cout << " parameter variables: (" << paramVarToString(var1) << ", " << paramVarToString(var2) << ") \nContinuing program..." << std::endl;
    }

  } else if (myParams.mission=="IceGiant"){
    //array of start values in units of years
    std::array<high_prec_t, 10> startTimes { 0., 10./9., 20./9., 30./9., 40./9., 50./9., 60./9., 70./9., 80./9., 90./9. };
    
    //iterate over start times. 10 40-day observation windows
    for (auto& stime: startTimes){
      stime *= Constants::yr; //convert to seconds
      status1 = gsl_integration_qag(&F1, stime, stime + 40.*Constants::day, epsilonAbsolute, epsilonRelative, 1000000, 6, workSpace, &result1, &error1);
      totalResult += 2*result1/sN(myParams.freqGW);
      totalError = sqrt(totalError*totalError + error1*error1);
    };


    //error handling
    if (status1==GSL_EROUND){
      std::cout << "Round off error detected. \n result = " << result1 << "  error = " << error1 << "  relative error1 = " << error1/result1 << std::endl;
      std::cout << " parameter variables: (" << paramVarToString(var1) << ", " << paramVarToString(var2) << ") \nContinuing program..." << std::endl;
    }

  } else {
    assert(("ERROR: Mission parameter not recognized", false));
  };

  //freeup memory used for workspace and return results
  gsl_integration_workspace_free(workSpace);
  return {totalResult, totalError};
}
*/


/* //array of start values in units of years
    std::array<high_prec_t, 10> startTimes { 0., 10./9., 20./9., 30./9., 40./9., 50./9., 60./9., 70./9., 80./9., 90./9. };

    high_prec_t Nsamples { 40.*Constants::day * myParams.freqGW * periodSamples };
    high_prec_t dx { (high_prec_t)40.*Constants::day/Nsamples};

    
    //iterate over start times. 10 40-day observation windows
    for (auto& stime: startTimes){
      
      stime *= Constants::yr; //convert to seconds
      //compute input arrays
      std::vector<high_prec_t> input1 { sampleFunction(function1, stime, stime + 40.*Constants::day, Nsamples, numThreads) };
      std::vector<high_prec_t> input2 { sampleFunction(function2, stime, stime + 40.*Constants::day, Nsamples, numThreads) };

      //caluclate ffts of input arrays
      std::vector<std::vector<high_prec_t>> fft1 { fourierTransform(input1, numThreads) };
      std::vector<std::vector<high_prec_t>> fft2 { fourierTransform(input2, numThreads) };

      //compute integrand by multiplying ffts together
      std::vector<high_prec_t> realPart{};
      std::vector<high_prec_t> complexPart{};
      for (int i{0}; i<Nsamples; ++i){
        realPart.push_back(fft1[i][0]*fft2[i][0] + fft1[i][1]*fft2[i][1]);
        complexPart.push_back(fft1[i][1]*fft2[i][0] - fft1[i][0]*fft2[i][1]);
      }

      //compute integral
      high_prec_t integral1 { myTrapSum(realPart, dx, numThreads) };
      high_prec_t integral1c { myTrapSum(complexPart, dx, numThreads) };
      if (integral1c>1e-10){
        std::cout << "Complex part of integral: " << integral1c << std::endl;
      }*/