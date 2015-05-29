#ifndef CONFIGPARAMSHOMOG_H
#define CONFIGPARAMSHOMOG_H

#include "ConfigParams.h"

namespace USACConfig
{
	// problem specific/data-related parameters: fundamental matrix
	struct Homog
	{
		Homog(): inputFilePath	      ("")			// leave blank if not using config file
		{}

		std::string			inputFilePath;
	};
}

class ConfigParamsHomog: public ConfigParams
{
public:
	// simple function to read in parameters from config file
	bool initParamsFromConfigFile(std::string& configFilePath);

	USACConfig::Homog homog;
};

#endif