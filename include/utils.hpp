#ifndef UTILS_H_INCLUDED
#define UTILS_H_INCLUDED

#include "Constants.hpp"
#include "VariablesParameters.hpp"
#include <iostream>
#include <chrono> // for std::chrono functions
#include <vector>
#include <string>
#include <functional>
using namespace std;


bool isClose(high_prec_t a, high_prec_t b, high_prec_t rel_tol = 1e-4, high_prec_t abs_tol = 0.0);

high_prec_t cosThetaS(high_prec_t,high_prec_t,high_prec_t);

template <typename T>
T thetaS(T thetaS, T phiS, T time){
  T costheta { cosThetaS(thetaS, phiS, time) };
  return std::acos(costheta);
}

template <typename T>
T thetaS(params& myParams, T time){
  T costheta { cosThetaS(myParams.thetaS, myParams.phiS, time) };
  return std::acos(costheta);
}

high_prec_t CosIota(params&);


template <typename T>
auto phiS(T thetaS, T phiS, T time){
  //azimuthal angle phi_S of the source in the detector frame
  T numerator { sqrt(3) * cos(thetaS) + sin(thetaS)*cos(Constants::omegaE*time - phiS) };
  T denominator { 2*sin(thetaS) * sin(Constants::omegaE * time - phiS) };

  return Constants::omegaE*time + atan2(numerator, denominator);
}

template <typename T>
auto phiS(params& myParams, T time){
  return phiS(myParams.thetaS, myParams.phiS, time);
}

high_prec_t psiS(high_prec_t, high_prec_t, high_prec_t, high_prec_t, high_prec_t, high_prec_t);
high_prec_t psiS(params&, high_prec_t);

high_prec_t Sqrt(high_prec_t);
high_prec_t Sin(high_prec_t);
high_prec_t Cos(high_prec_t);
high_prec_t Power(high_prec_t, high_prec_t);
high_prec_t ArcTan(high_prec_t,high_prec_t);
high_prec_t ArcTan(high_prec_t);

std::string paramVarToString(ParameterVariables var);

high_prec_t cassiniPSD(params&, high_prec_t);

//timer class taken from --> https://www.learncpp.com/cpp-tutorial/timing-your-code/
class Timer
{
private:
	// Type aliases to make accessing nested type easier
	using Clock = std::chrono::steady_clock;
	using Second = std::chrono::duration<high_prec_t, std::ratio<1> >;

	std::chrono::time_point<Clock> m_beg { Clock::now() };

public:
	void reset()
	{
		m_beg = Clock::now();
	}

	high_prec_t elapsed() const
	{
		return std::chrono::duration_cast<Second>(Clock::now() - m_beg).count();
	}
};



#ifdef COMPILE_AS_PROGRAM
#include <eigen3/Eigen/Dense> //for matricies and related methods

using Eigen::MatrixXd;
class FisherMatrix{
private:
	MatrixXd m_matrix;
	MatrixXd m_errors;

	bool m_inverseCalculated{false};
	MatrixXd m_inverse;

public:
	unsigned int Rows;
	unsigned int Cols;

	//default constructor
	FisherMatrix() = default;

	//FisherMatrix constructor
	FisherMatrix(unsigned int x): m_matrix(x,x), m_errors(x,x), m_inverse(x,x){Rows = x; Cols = x;};

	//FisherMatrix constructor
	FisherMatrix(MatrixXd matrix) {m_matrix << matrix;};

	//FisherMatrix list constructor
	//FisherMatrix(high_prec_t inputMatrix): m_matrix(inputMatrix){};

	//overloaded operator[] to get specific row
	//std::vector<high_prec_t>& operator[](unsigned int x);

	//methods to access matrix, errors, and inverse
	MatrixXd matrix() const {return m_matrix;};
	MatrixXd errors() const {return m_errors;};
	MatrixXd inverse() const {
		assert(m_inverseCalculated);
		return m_inverse;
	}

	void printMatrix() const { std::cout << m_matrix << "\n";};
	void printErrors() const { std::cout << m_errors << "\n";};

	//add to errors
	void addError(double value, std::vector<unsigned int> location){m_errors(location[0], location[1])=value;};

	//Calculate inverse
	void calculateInverse(){ m_inverse = m_matrix.inverse(); m_inverseCalculated=true;};

	//overload operator() to access matrix elements from class
	double& operator()(unsigned x, unsigned y);

	//overload multiplication operator
	void operator*(double val);
	
	//print inverse if calculated. if not, first calculate, then print
	void printInverse();

	std::vector<std::vector<double>> matrixRaw();

	std::vector<std::vector<double>> errorsRaw();

	std::vector<std::vector<double>> inverseRaw();
};


std::vector<high_prec_t> sampleFunction(std::function<high_prec_t(high_prec_t)>, high_prec_t, high_prec_t, int, int);
high_prec_t myTrapSum(std::vector<high_prec_t>&, high_prec_t, int);
high_prec_t myTrapIntegral(std::function<high_prec_t(high_prec_t)>, high_prec_t, high_prec_t, int, int);
std::vector<std::vector<high_prec_t>> fourierTransform(std::vector<high_prec_t>, int);

#endif //COMPILE_AS_PROGRAM


#endif //UTILS_H_INCLUDED
