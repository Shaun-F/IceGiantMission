#include "utils.hpp"
#include "VariablesParameters.hpp"
#include <omp.h>
#include <cmath>
#include <string>
#include <numeric>
#include <complex.h>
#include <fftw3.h>
using namespace std;


bool isClose(high_prec_t a, high_prec_t b, high_prec_t rel_tol, high_prec_t abs_tol){
  high_prec_t threshold = max( rel_tol * max(abs(a), abs(b)), abs_tol);
  return abs(a - b) <= threshold;
}

high_prec_t cosThetaS(high_prec_t thetaS, high_prec_t phiS, high_prec_t time){
  // cos(theta_S) in the detector frame
  high_prec_t val {sqrt(3)/2};
  high_prec_t term1 { 0.5*cos(thetaS) };
  high_prec_t term2 { val * sin(thetaS) * cos(Constants::omegaE * time - phiS) };

  return term1 - term2;
}

high_prec_t CosIota(params& params){
  return  cos(params.thetaL) * cos(params.thetaS) + sin(params.thetaL)*sin(params.thetaS)*cos(params.phiL - params.phiS);
}

high_prec_t psiS(high_prec_t thetaL, high_prec_t thetaS, high_prec_t phiL, high_prec_t phiS, high_prec_t cosIota, high_prec_t time){
  //polarization angle psi_S of wavefront in detector frame
  high_prec_t val {sqrt(3)/2};
  high_prec_t Lz { 0.5 * cos(thetaL) - val * sin(thetaL) * cos(Constants::omegaE * time - phiL) };

  high_prec_t nLxz_term1 { 0.5 * sin(thetaL) * sin(thetaS) * sin(phiL - phiS) };
  high_prec_t nLxz_term2 { val * cos(Constants::omegaE*time) * ( cos(thetaL) * sin(thetaS)*sin(phiS) - cos(thetaS)*sin(thetaL)*sin(phiL)) };
  high_prec_t nLxz_term3 { val * sin(Constants::omegaE*time) * ( cos(thetaS) * sin(thetaL)*cos(phiL) - cos(thetaL)*sin(thetaS)*cos(phiS)) };
  high_prec_t nLxz = { nLxz_term1 - nLxz_term2 - nLxz_term3 };

  high_prec_t numerator { Lz - cosIota*cosThetaS(thetaS, phiS, time) };
  return atan2(numerator, nLxz);
}

high_prec_t psiS(params& params, high_prec_t time){
  high_prec_t thetaL{params.thetaL};
  high_prec_t thetaS{params.thetaS};
  high_prec_t phiL{params.phiL};
  high_prec_t phiS{params.phiS};
  high_prec_t cosIota{CosIota(params)};
  high_prec_t result = {};
  high_prec_t numerator {};
  if (params.mission=="LISA"){
      //polarization angle psi_S of wavefront in detector frame for LISA
      high_prec_t val {sqrt(3)/2};
      high_prec_t Lz { 0.5 * cos(thetaL) - val * sin(thetaL) * cos(Constants::omegaE * time - phiL) };

      high_prec_t nLxz_term1 { 0.5 * sin(thetaL) * sin(thetaS) * sin(phiL - phiS) };
      high_prec_t nLxz_term2 { val * cos(Constants::omegaE*time) * ( cos(thetaL) * sin(thetaS)*sin(phiS) - cos(thetaS)*sin(thetaL)*sin(phiL)) };
      high_prec_t nLxz_term3 { val * sin(Constants::omegaE*time) * ( cos(thetaS) * sin(thetaL)*cos(phiL) - cos(thetaL)*sin(thetaS)*cos(phiS)) };
      high_prec_t nLxz = { nLxz_term1 - nLxz_term2 - nLxz_term3 };
      high_prec_t numerator {Lz - cosIota*cosThetaS(thetaS, phiS, time)};
      result = atan2(numerator, nLxz);
  } else if (params.mission =="IceGiant"){
      //polarization angle psi_S of wavefront in detector frame for Ice Giant mission
      high_prec_t val {sqrt(3)/2};
      high_prec_t Lz { 0.5 * cos(thetaL) - val * sin(thetaL) * cos(params.ig_direction- phiL) };

      high_prec_t nLxz_term1 { 0.5 * sin(thetaL) * sin(thetaS) * sin(phiL - phiS) };
      high_prec_t nLxz_term2 { val * cos(params.ig_direction) * ( cos(thetaL) * sin(thetaS)*sin(phiS) - cos(thetaS)*sin(thetaL)*sin(phiL)) };
      high_prec_t nLxz_term3 { val * sin(params.ig_direction) * ( cos(thetaS) * sin(thetaL)*cos(phiL) - cos(thetaL)*sin(thetaS)*cos(phiS)) };
      high_prec_t nLxz = { nLxz_term1 - nLxz_term2 - nLxz_term3 };
      high_prec_t numerator { Lz - cosIota*cos(params.thetaS)};

      result = atan2(numerator, nLxz);
  } else {
      assert(("Mission parameter not known. ", false));
  };

  return result;
}

high_prec_t Sqrt(high_prec_t x){
  return sqrt(x);
}
high_prec_t Sin(high_prec_t x){
  return sin(x);
}
high_prec_t Cos(high_prec_t x){
  return cos(x);
}
high_prec_t Power(high_prec_t x, high_prec_t y){
  return pow(x,y);
}
high_prec_t ArcTan(high_prec_t x, high_prec_t y){
  return atan2(y,x);
}
high_prec_t ArcTan(high_prec_t x){
  return atan(x);
}

high_prec_t cassiniPSD(params& myParams, high_prec_t freq){
  high_prec_t result;
  high_prec_t PSDlevel { myParams.PSDlevel };
  high_prec_t alpha1 { -2 };
  high_prec_t alpha2 { -1./2. };
  high_prec_t alpha3 { 2 };
  high_prec_t f1 { 7.07e-5 };
  high_prec_t f2 { 2.5e-5 };
  high_prec_t f3 { 7.95e-1 };

  result = PSDlevel * ( std::pow(freq/f1, alpha1) + std::pow(freq/f2, alpha2) + std::pow(freq/f3, alpha3) );

  return result;
};

std::string paramVarToString(ParameterVariables var){
  std::string output{};
  switch(var){
    case ParameterVariables::K:
      output="K";
      break;
    case ParameterVariables::P:
      output="P";
      break;
    case ParameterVariables::phiP:
      output="phiP";
      break;
    case ParameterVariables::thetaS:
      output="thetaS";
      break;
    case ParameterVariables::phiS:
      output="phiS";
      break;
    case ParameterVariables::thetaL:
      output="thetaL";
      break;
    case ParameterVariables::phiL:
      output="phiL";
      break;
    case ParameterVariables::lnA:
      output="lnA";
      break;
    case ParameterVariables::f1:
      output="f1";
      break;
    case ParameterVariables::f0:
      output="f0";
      break;
  }

  return output;
};


#ifdef COMPILE_AS_PROGRAM

high_prec_t myTrapSum(std::vector<high_prec_t>& arr, high_prec_t dx, int numThreads=1){
    
    high_prec_t integrationResult { 2*std::reduce(++arr.cbegin(), --arr.cend()) };
    integrationResult += arr[0] + arr[arr.size()-1];
    integrationResult *= (1./2.)*dx;

    return integrationResult;
}

std::vector<high_prec_t> sampleFunction(std::function<high_prec_t(high_prec_t)> func, high_prec_t start, high_prec_t stop, int nPoints, int numThreads=1){
  high_prec_t dx { (high_prec_t)(stop-start) / nPoints };
  std::vector<high_prec_t> funcSample;
  std::vector<std::vector<high_prec_t>> funcSamplePrivate{};
  funcSamplePrivate.resize(numThreads);

  #pragma omp parallel num_threads(numThreads)
	{
		#pragma omp for nowait
		for (int i=0; i<nPoints; ++i){
			int threadID {omp_get_thread_num()};
			high_prec_t abscissa { start + i*dx };
			funcSamplePrivate[threadID].push_back(func(abscissa));			
		};
	}
	for (int i{0}; i<numThreads; ++i){
		funcSample.insert(funcSample.end(), funcSamplePrivate[i].begin(), funcSamplePrivate[i].end());
	}

  return funcSample;
}

high_prec_t myTrapIntegral(std::function<high_prec_t(high_prec_t)> func, high_prec_t start, high_prec_t stop, int nPoints, int numThreads=1){
  high_prec_t dx { (high_prec_t)(stop-start) / nPoints };
  std::vector<high_prec_t> funcSample { sampleFunction(func, start, stop, nPoints, numThreads) };

	high_prec_t integrationResult { myTrapSum(funcSample, dx, numThreads) };

	return integrationResult;
}

//computes fourier transform of REAL array
std::vector<std::vector<high_prec_t>> fourierTransform(std::vector<high_prec_t> input, int num_threads=1){
  

  //setup initial containers
  std::vector<std::vector<high_prec_t>> output {};
  int length { (int)input.size() };
  fftw_complex *in, *out;
  
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * length);
  
  //populate input array
  #pragma omp parallel for num_threads(num_threads)
  for (int i=0; i<length; ++i){
    in[i][0] = input[i];
    in[i][1] = 0.;
  }
  
  //build fft plan
  fftw_plan fftPlan { fftw_plan_dft_1d(length, in, out, FFTW_FORWARD, FFTW_ESTIMATE) };
  
  //execute fft computation
  fftw_execute(fftPlan);
	
  //build return container
  for (int i{0}; i<length; ++i){
    output.push_back( 
          { 
            (double)out[i][0]/std::sqrt(length), //note normalization factor. such that fourier inverse of fourier transform is exactly original array.
						(double)out[i][1]/std::sqrt(length) 
					  }
					);
  }

  fftw_destroy_plan(fftPlan);

  //return fft
  return output;
}



double& FisherMatrix::operator()(unsigned x, unsigned y){ 
  return m_matrix(x,y);
};

void FisherMatrix::operator*(double val){
  m_matrix *= val; 
  if (m_inverseCalculated){
    calculateInverse();
  }
};

void FisherMatrix::printInverse(){
	  if (!m_inverseCalculated){
	    calculateInverse();
	  };
	  std::cout << m_inverse << "\n";
	};

std::vector<std::vector<double>> FisherMatrix::matrixRaw(){
		std::vector<std::vector<double>> rawData;
		rawData.resize(Rows);
		for (auto& el: rawData){
			el.resize(Cols);
		};

		for (int row{0}; row<Rows; ++row){
			for(int col{0}; col<Cols; ++col){
				rawData[row][col] = m_matrix(row,col);
			};
		};
		return rawData;
	};

std::vector<std::vector<double>> FisherMatrix::errorsRaw(){
		std::vector<std::vector<double>> rawData;
		rawData.resize(Rows);
		for (auto& el: rawData){
			el.resize(Cols);
		};

		for (int row{0}; row<Rows; ++row){
			for(int col{0}; col<Cols; ++col){
				rawData[row][col] = m_errors(row,col);
			};
		};
		return rawData;
	};

std::vector<std::vector<double>> FisherMatrix::inverseRaw(){
		assert(m_inverseCalculated);
		std::vector<std::vector<double>> rawData;
		rawData.resize(Rows);
		for (auto& el: rawData){
			el.resize(Cols);
		};

		for (int row{0}; row<Rows; ++row){
			for(int col{0}; col<Cols; ++col){
				rawData[row][col] = m_inverse(row,col);
			};
		};
		return rawData;
	};

#endif //COMPILE_AS_PROGRAM

//constructor
/*
FisherMatrix::FisherMatrix(unsigned int x){
  m_matrix.resize(x,x);
  m_errors.resize(x,x);
  m_inverse.resize(x,x);

};
*/

//overload operator[]
/*
std::vector<high_prec_t>& FisherMatrix::operator[](unsigned int x){
  return m_matrix.at(x);
}
*/
