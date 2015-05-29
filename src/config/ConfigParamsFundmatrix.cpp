#include "ConfigParamsFundmatrix.h"
#include "ConfigFileReader.h"

bool ConfigParamsFund::initParamsFromConfigFile(std::string &configFilePath)
{
	// initialize the config file reader with the input config file
	ConfigFileReader cfr;
	try 
	{
		cfr.readFile(configFilePath.c_str());
	}
	catch(const libconfig::FileIOException)
	{
		std::cerr << "I/O error while reading file." << std::endl;
		return false;
	}
	catch(const libconfig::ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
			<< " - " << pex.getError() << std::endl;
		return false;
	}

	// now read in parameters 
	try
	{
		// usac parameters
		if (!initUSACParamsFromConfigFile(cfr)) 
		{
			
			std::cerr << "Error reading USAC parameters from config file." << std::endl;
			return false;
		}

		// fundmatrix parameters
		if ( !cfr.getTopLevelParameterValue("problem_specific", "input_file_path", fund.inputFilePath) )
		{
			return false;
		}

		if ( common.testDegeneracy )
		{
			if ( !(cfr.getTopLevelParameterValue("problem_specific", "h_degen_threshold", fund.hDegenThreshold) &&
				   cfr.getTopLevelParameterValue("problem_specific", "max_upgrade_samples", fund.maxUpgradeSamples) ))
			{
				return false;
			}
		}

		std::string decomp_type;
		if ( !cfr.getTopLevelParameterValue("problem_specific", "matrix_decomposition", decomp_type) )
		{
			return false;
		}
		if (decomp_type.compare("LU") == 0)
		{
			fund.decompositionAlg = USACConfig::DECOMP_LU;
		}
		else if (decomp_type.compare("QR") == 0)
		{
			fund.decompositionAlg = USACConfig::DECOMP_QR;
		}
		else
		{
			std::cerr << "Matrix decomposition " << decomp_type << " not supported" << std::endl;
			return false;
		}
	}
	catch(...)
	{
		return false;
	}

	return true;
}