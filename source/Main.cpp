#include "utils.hpp"
#include "LISA.hpp"
#include "Binary.hpp"
#include "Constants.hpp"
#include "VariablesParameters.hpp"
#include "json.h" //contains json type
#include <omp.h> //parallelism
#include <fstream> //for data outputting
#include <iostream> //outputting to terminal
#include <filesystem> //for getting directory info

//for parsing ini file
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


//setup namespace for json class
using nlohmann::json;

//forward declarations for functions used in main()
void executeAnalysis(params&, Binary&, json&);
void combineInformation(params&, analysisParams&, Binary&, json&);
void saveJson(json&, std::ostream&);
json readJson(std::ifstream&);
void getParameters(params&, const std::string_view);
void getAnalysisParameters(analysisParams&, const std::string_view);
void printParameters(const params&);
void printAnalysisParameters(const analysisParams&);
template <typename T> std::vector<T> to_array(const std::string& s);


int main(int argc, char** argv){
  /*
  NOTES:
    there seems to be a sign error in the GW polarization amplitudes in ig_binary (see self.gw_amplitude in ig_binary constructor)
	we use T=15000 for the two way light travel time, while the default in ig_binary is T=30000
  */

	//reading in parameters from file
	params myParams;
	getParameters(myParams, argv[1]);
 	
	//print parameters to screen
	printParameters(myParams);

	//read analysis-specific parameters
	analysisParams AnalysisParameters;
	getAnalysisParameters(AnalysisParameters, argv[1]);

	//print analysis parameters to screen
	printAnalysisParameters(AnalysisParameters);

	//output current precision of code to screen
	std::cout << "Using precision = " << std::numeric_limits<high_prec_t>::digits10 << std::endl;

	//get current directory
	auto dataPath {std::filesystem::current_path()}; //get current working directory

	//set file output
	std::string modePrefix { myParams.mode };
	std::string dataFileName { modePrefix + "_CBPJointSearchData.json" };
	dataPath /= dataFileName;
	std::cout << "Data will be stored in the file: " << dataPath << std::endl;
	
	//instantiate json object for data saving
	std::ifstream dataFile {dataPath};
	json jobj = readJson(dataFile);

	//timer object
	Timer globalTimer{};

	//iterate over periods for CBP
	for (high_prec_t& period: AnalysisParameters.PeriodList){
		
		//set new period
		myParams.P = period*Constants::yr;

		std::cout << "Running analysis with CBP period of "<< myParams.P/Constants::yr << " years. " << std::endl;

		//instantiate binary object
		Binary myBinary{myParams};

		//run analysis for LISA mission
		myParams.mission = "LISA";
		executeAnalysis(myParams, myBinary, jobj);

		//run analysis for IceGiant mission
		myParams.mission = "IceGiant";
		executeAnalysis(myParams, myBinary, jobj);

		//combine information matrices to get joint-search uncertainties
		combineInformation(myParams, AnalysisParameters, myBinary, jobj);

		//save data to file after each iteration
		std::ofstream DataFile{dataPath};
		saveJson(jobj, DataFile);

	};//endfor
	
	std::cout << "Program complete. Time spent on iterations: " << globalTimer.elapsed() << std::endl;

  return 0;
}

void executeAnalysis(params& myParams, Binary& myBinary, json& jobj){

	//instantiate timer object to track execution times
	Timer myTimer{};

	unsigned int FisherMatrixSize{};
	if (myParams.mode=='s'){
		FisherMatrixSize=3;
	} else if (myParams.mode=='m'){
		FisherMatrixSize=9;
	} else if (myParams.mode=='l'){
		FisherMatrixSize=static_cast<int>(ParameterVariables::NUM_Variables);
	} else {
		assert(("ERROR: Mode parameter not recognized", false));
	};

	//instantiate fisher matrix object to save matrix to.
	FisherMatrix myFisherMatrix{FisherMatrixSize};

	//reset timer to 0
	myTimer.reset();

	//generate the fisher matrix
	std::cout << "Generating fisher matrix for "<<myParams.mission<<" mission..."<<std::endl;
	high_prec_t LISA_SNR;
	myBinary.genFisherMatrix(myParams, myFisherMatrix,LISA_SNR, (int)omp_get_max_threads());
	myFisherMatrix.calculateInverse();
	std::cout << "Fisher matrix generated in " << myTimer.elapsed() << " seconds." << std::endl;

	//save data to json object

	//calculate raw nested vectors from Eigen::MatrixXd objects
	std::vector<std::vector<double>> myMatrixData{myFisherMatrix.matrixRaw()};
	myFisherMatrix.calculateInverse();
	std::vector<std::vector<double>> myInverseData{myFisherMatrix.inverseRaw()};

	//formulate relative uncertainties
	//get array of parameter values in order of ParameterVariables struct
	std::vector<double> myParameters {
									(double)myBinary.FrequencyK(myParams), 
									(double)myParams.P, 
									(double)myParams.phiP, 
									(double)myParams.thetaS, 
									(double)myParams.phiS, 
									(double)myParams.thetaL, 
									(double)myParams.phiL, 
									(double)myParams.a0, 
									(double)myParams.f1, 
									(double)myParams.freqGW
									};

	//vector of relative uncertainties
	std::vector<double> relUnc{};
	std::vector<double> CovarianceMatrixDiag{};
	for (int i{0}; i<myFisherMatrix.Rows; ++i){
		CovarianceMatrixDiag.push_back(myInverseData[i][i]);
	}
	for (int i{0}; i<myFisherMatrix.Rows; ++i){
		relUnc.push_back(sqrt(CovarianceMatrixDiag[i])/myParameters[i]);
	};


	jobj[std::to_string(myParams.P)][myParams.mission]["SNR"] = LISA_SNR;
	jobj[std::to_string(myParams.P)][myParams.mission]["fisherMatrix"] = myMatrixData;
	jobj[std::to_string(myParams.P)][myParams.mission]["varianceMatrix"] = myInverseData;
	jobj[std::to_string(myParams.P)][myParams.mission]["relativeUncertainties"] = relUnc;
};

void combineInformation(params& myParams, analysisParams& AnalysisParams, Binary& myBinary, json& jobj){

	std::vector<std::vector<double>> LISAfisher { jobj[std::to_string(myParams.P)]["LISA"]["fisherMatrix"].template get<std::vector<std::vector<double>>>() }; //convert from json object types to std::vector
	unsigned int FisherMatrixSize { static_cast<unsigned int>(LISAfisher.size()) };

	//calculate total SNR
	double LISASNR {jobj[std::to_string(myParams.P)]["LISA"]["SNR"].template get<double>() };
	double totalSNR = sqrt(1+myParams.relativeSNR*myParams.relativeSNR) * LISASNR;

	FisherMatrix TotalMatrix{ FisherMatrixSize };

	//add data to Total Fisher Matrix
	for (unsigned int row{0}; row<jobj[std::to_string(myParams.P)]["LISA"]["fisherMatrix"].size(); ++row){
		for (unsigned int col{0}; col<jobj[std::to_string(myParams.P)]["LISA"]["fisherMatrix"][row].size(); ++col){
			double lisaElement { jobj[std::to_string(myParams.P)]["LISA"]["fisherMatrix"][row][col] };
			double icegiantElement { jobj[std::to_string(myParams.P)]["IceGiant"]["fisherMatrix"][row][col] };
			TotalMatrix(row,col) = lisaElement + icegiantElement; //add fisher matrices elementwise
		}
	};

	//calculate variance matrix
	TotalMatrix.calculateInverse();
	std::vector<std::vector<double>> myInverseData{TotalMatrix.inverseRaw()};
	std::vector<double> TotalVariances{};
	for (int i{0}; i<TotalMatrix.Rows; ++i){
		TotalVariances.push_back(myInverseData[i][i]);
	}

	//add relative uncertainties
	std::vector<double> myParameters {
									(double)myBinary.FrequencyK(myParams), 
									(double)myParams.P, 
									(double)myParams.phiP, 
									(double)myParams.thetaS, 
									(double)myParams.phiS, 
									(double)myParams.thetaL, 
									(double)myParams.phiL, 
									(double)myParams.a0, 
									(double)myParams.f1, 
									(double)myParams.freqGW
									};	
	std::vector<double> relUnc{};
	for (int i{0}; i<TotalMatrix.Rows; ++i){
		relUnc.push_back(sqrt(TotalVariances[i])/myParameters[i]);
	};

	jobj[std::to_string(myParams.P)]["Total"]["fisherMatrix"] = TotalMatrix.matrixRaw();
	jobj[std::to_string(myParams.P)]["Total"]["relativeUncertainties"] = relUnc;
	jobj[std::to_string(myParams.P)]["Total"]["varianceMatrix"] = TotalVariances;
	jobj[std::to_string(myParams.P)]["Total"]["SNR"] = totalSNR;
};

void saveJson(json& jobj, std::ostream& outFile){
	outFile << std::setw(4)<<jobj << std::endl;
}

json readJson(std::ifstream& readFile){
	bool fileExists { readFile.is_open() };
	
	if (!fileExists){
		json jobj; //create empty json object
		return jobj;
	} else {
		json jobj = json::parse(readFile);
		return jobj;
	};
}

void getParameters(params& out, const std::string_view fileName){
	//get current directory
	auto dataPath {std::filesystem::current_path()}; //get current working directory

	//set file output
	std::string dataFileName { fileName };
	dataPath /= dataFileName;
	
	boost::property_tree::ptree pt; //property tree object
	boost::property_tree::ini_parser::read_ini(dataPath, pt); //parse the ini file

	out.mission			= pt.get("parameters.mission", "LISA");
	out.thetaL 			= pt.get("parameters.thetaL", 0.);
	out.phiL			= pt.get("parameters.phiL", 0.);
	out.thetaS 			= pt.get("parameters.thetaS", 0.);
	out.phiS			= pt.get("parameters.phiS", 0.);
	out.LISAAlpha		= pt.get("parameters.LISAAlpha", 2);	
	out.M1				= pt.get("parameters.M1", 0.);
	out.M2				= pt.get("parameters.M2", 0.);
	out.MP				= pt.get("parameters.MP", 0.);
	out.Tobs			= pt.get("parameters.Tobs", 4.0*Constants::yr);
	out.Larm			= pt.get("parameters.Larm", 2.5e9);
	out.NC				= pt.get("parameters.NC", 2);
	out.thetaP			= pt.get("parameters.thetaP", 0.);
	out.phiP			= pt.get("parameters.phiP", 0.);
	out.mode			= pt.get("parameters.mode", "s")[1];
	out.freqGW			= pt.get("parameters.freqGW", 0.);
	out.sourceDistance	= pt.get("parameters.sourceDistance", 0.);
	out.DerivativeDelta	= pt.get("parameters.DerivativeDelta", 1e-6);
	out.ig_direction	= pt.get("parameters.ig_direction", 0.);
	out.lightTwoWayTime	= pt.get("parameters.lightTwoWayTime", 15000);
	out.periodSamples	= pt.get("parameters.periodSamples", 35);	
	out.PSDlevel		= pt.get("parameters.PSDlevel", 9e-26);
	out.allanDeviation	= pt.get("parameters.allanDeviation", 3e-15);
	out.relativeSNR		= pt.get("parameters.relativeSNR", 1e-1);

}

void getAnalysisParameters(analysisParams& out, const std::string_view fileName){
	//get current directory
	auto dataPath {std::filesystem::current_path()}; //get current working directory

	//set file output
	std::string dataFileName { fileName };
	dataPath /= dataFileName;
	
	boost::property_tree::ptree pt; //property tree object
	boost::property_tree::ini_parser::read_ini(dataPath, pt); //parse the ini file

	out.PeriodList		= to_array<high_prec_t>(pt.get<std::string>("analysis.PeriodList"));
}

void printParameters(const params& myParams){
	std::string header {"#############################################\n"};
	std::string footer { header };

	std::string out { 	"thetaL: " + std::to_string(myParams.thetaL) +
						"\nphiL: " + std::to_string(myParams.phiL) + 
						"\nthetaS: " + std::to_string(myParams.thetaS) + 
						"\nphiS: " + std::to_string(myParams.phiS) + 
						"\nlisaAlpha: " + std::to_string(myParams.LISAAlpha) + 
						"\nM1: " + std::to_string(myParams.M1) + 
						"\nM2: " + std::to_string(myParams.M2) + 
						"\nMP: " + std::to_string(myParams.MP) + 
						"\nTobs: " + std::to_string(myParams.Tobs) + 
						"\nLarm: " + std::to_string(myParams.Larm) + 
						"\nNC: " + std::to_string(myParams.NC) + 
						"\nthetaP: " + std::to_string(myParams.thetaP) + 
						"\nphiP: " + std::to_string(myParams.phiP) + 
						"\nmode: " + myParams.mode + 
						"\nfreqGW: " + std::to_string(myParams.freqGW) + 
						"\nsourceDistance: " + std::to_string(myParams.sourceDistance) + 
						"\nderivativeDelta: " + std::to_string(myParams.DerivativeDelta) + 
						"\nig_direction: " + std::to_string(myParams.ig_direction) + 
						"\nlightTwoWayTime: " + std::to_string(myParams.lightTwoWayTime) +
						"\nperiodSamples: " + std::to_string(myParams.periodSamples) 	+
						"\nrelativeSNR: " + std::to_string(myParams.relativeSNR)
					};

	std::cout << header + header + header << out << "\n" <<footer + footer + footer << std::endl;
}

void printAnalysisParameters(const analysisParams& anPars){
	std::string header {"#############################################\n"};
	std::string footer { header };
	std::string out {"Periodlist: {"};
	for (auto el: anPars.PeriodList){
		out += std::to_string(el) + ", ";
	}

	std::cout << header + header + header << out << "\n" <<footer + footer + footer << std::endl;
}

template<typename T>
std::vector<T> to_array(const std::string& s)
{
  std::vector<T> result;
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, ',')) result.push_back(boost::lexical_cast<T>(item));
  return result;
}
