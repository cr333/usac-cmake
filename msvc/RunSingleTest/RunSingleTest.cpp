#define NOMINMAX
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <windows.h>
#include "config\ConfigParams.h"
#include "estimators\FundMatrixEstimator.h"
#include "estimators\HomogEstimator.h"

// helper functions
bool readCorrsFromFile(std::string& inputFilePath, std::vector<double>& pointData, unsigned int& numPts)
{
	// read data from from file
	std::ifstream infile(inputFilePath.c_str());
	if (!infile.is_open())
	{
		std::cerr << "Error opening input points file: " << inputFilePath << std::endl;
		return false;
	}
	infile >> numPts;
	pointData.resize(6*numPts);
	for (unsigned int i = 0; i < numPts; ++i)
	{
		infile >> pointData[6*i] >> pointData[6*i+1] >> pointData[6*i+3] >> pointData[6*i+4];
		pointData[6*i+2] = 1.0;
		pointData[6*i+5] = 1.0;
	}
	infile.close();
	return true;
}

bool readPROSACDataFromFile(std::string& sortedPointsFile, unsigned int numPts, std::vector<unsigned int>& prosacData)
{
	std::ifstream infile(sortedPointsFile.c_str());
	if (!infile.is_open())
	{
		std::cerr << "Error opening sorted indices file: " << sortedPointsFile << std::endl;
		return false;
	}
	for (unsigned int i = 0; i < numPts; ++i)
	{
		infile >> prosacData[i];
	}
	infile.close();
	return true;
}


// ---------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv)
{
	// check command line args
	if (argc < 3)
	{
		std::cerr << "Usage: RunSingleTest <estimation problem> <config file>" << std::endl;
		std::cerr << "\t<estimation problem>: 0 (fundamental matrix), 1 (homography)" << std::endl;
		std::cerr << "\t<config file>: full path to configuration file" << std::endl;
		return(EXIT_FAILURE);
	}
	int estimation_problem = atoi(argv[1]);
	std::string cfg_file_path = argv[2];

	// seed random number generator
	srand((unsigned int)time(NULL));

	// initialize the appropriate robust estimation problem
	if (estimation_problem == 0)
	{
		// ------------------------------------------------------------------------
		// initialize the fundamental matrix estimation problem
		ConfigParamsFund cfg;
		if ( !cfg.initParamsFromConfigFile(std::string(cfg_file_path)) )
		{
			std::cerr << "Error during initialization" << std::endl;
			return(EXIT_FAILURE);
		}
		FundMatrixEstimator* fund = new FundMatrixEstimator;
		fund->initParamsUSAC(cfg);

		// read in input data points
		std::vector<double> point_data;
		if ( !readCorrsFromFile(cfg.fund.inputFilePath, point_data, cfg.common.numDataPoints) )
		{
			return(EXIT_FAILURE);
		}

		// read in prosac data if required
		std::vector<unsigned int> prosac_data;
		if (cfg.common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
		{
			prosac_data.resize(cfg.common.numDataPoints);
			if ( !readPROSACDataFromFile(cfg.prosac.sortedPointsFile, cfg.common.numDataPoints, prosac_data) )
			{
				return(EXIT_FAILURE);
			}
			cfg.prosac.sortedPointIndices = &prosac_data[0];
		} else {
			cfg.prosac.sortedPointIndices = NULL;
		}

		// set up the fundamental matrix estimation problem

		fund->initDataUSAC(cfg);
		fund->initProblem(cfg, &point_data[0]);
		if (!fund->solve())
		{
			return(EXIT_FAILURE);
		}

		// write out results
		size_t pos = (cfg.fund.inputFilePath).find_last_of("/\\");
		std::string working_dir = (cfg.fund.inputFilePath).substr(0, pos);
		std::ofstream outmodel((working_dir + "\\F.txt").c_str());
		for (unsigned int i = 0; i < 3; ++i)
		{
			for (unsigned int j = 0; j < 3; ++j)
			{
				outmodel << fund->final_model_params_[3*i+j] << " ";
			}
		}
		outmodel.close();
		std::ofstream outinliers((working_dir + "\\inliers.txt").c_str());
		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
		{
			outinliers << fund->usac_results_.inlier_flags_[i] << std::endl;
		}
		outinliers.close();

		// clean up
		point_data.clear();
		prosac_data.clear();
		fund->cleanupProblem();
		delete fund;

	} else if (estimation_problem == 1) {
		// ------------------------------------------------------------------------
		// initialize the homography estimation problem
		ConfigParamsHomog cfg;
		if ( !cfg.initParamsFromConfigFile(std::string(cfg_file_path)) )
		{
			std::cerr << "Error during initialization" << std::endl;
			return(EXIT_FAILURE);
		}

		HomogEstimator* homog = new HomogEstimator;
		homog->initParamsUSAC(cfg);

		// read in input data points
		std::vector<double> point_data;
		if ( !readCorrsFromFile(cfg.homog.inputFilePath, point_data, cfg.common.numDataPoints) )
		{
			return(EXIT_FAILURE);
		}

		// read in prosac data if required
		std::vector<unsigned int> prosac_data;
		if (cfg.common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
		{
			prosac_data.resize(cfg.common.numDataPoints);
			if ( !readPROSACDataFromFile(cfg.prosac.sortedPointsFile, cfg.common.numDataPoints, prosac_data) )
			{
				return(EXIT_FAILURE);
			}
			cfg.prosac.sortedPointIndices = &prosac_data[0];
		} else {
			cfg.prosac.sortedPointIndices = NULL;
		}

		// set up the homography estimation problem
		homog->initDataUSAC(cfg);
		homog->initProblem(cfg, &point_data[0]);
		if (!homog->solve())
		{
			return(EXIT_FAILURE);
		}

		// write out results
		size_t pos = (cfg.homog.inputFilePath).find_last_of("/\\");
		std::string working_dir = (cfg.homog.inputFilePath).substr(0, pos);
		std::ofstream outmodel((working_dir + "\\H.txt").c_str());
		for (unsigned int i = 0; i < 3; ++i)
		{
			for (unsigned int j = 0; j < 3; ++j)
			{
				outmodel << homog->final_model_params_[3*i+j] << " ";
			}
		}
		outmodel.close();
		std::ofstream outinliers((working_dir + "\\inliers.txt").c_str());
		for (unsigned int i = 0; i < cfg.common.numDataPoints; ++i)
		{
			outinliers << homog->usac_results_.inlier_flags_[i] << std::endl;
		}
		outinliers.close();

		// clean up
		point_data.clear();
		prosac_data.clear();
		homog->cleanupProblem();
		delete homog;
	} else {
		std::cout << "Estimation problem currently not implemented" << std::endl;
	}

	return(EXIT_SUCCESS);
}

