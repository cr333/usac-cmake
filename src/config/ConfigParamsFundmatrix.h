#ifndef CONFIGPARAMSFUND_H
#define CONFIGPARAMSFUND_H

#include "ConfigParams.h"

namespace USACConfig
{
	enum MatrixDecomposition      {DECOMP_QR, DECOMP_LU};

	// problem specific/data-related parameters: fundamental matrix
	struct Fund
	{
		Fund(): decompositionAlg  (DECOMP_QR),
			hDegenThreshold       (0.0),
			maxUpgradeSamples     (500),
			inputFilePath	      ("")			// leave blank if not using config file
		{}

		MatrixDecomposition decompositionAlg;
		double				hDegenThreshold;
		unsigned int		maxUpgradeSamples;
		std::string			inputFilePath;
	};
}

class ConfigParamsFund: public ConfigParams
{
public:
	// simple function to read in parameters from config file
	bool initParamsFromConfigFile(std::string& configFilePath);

	USACConfig::Fund fund;
};

#endif