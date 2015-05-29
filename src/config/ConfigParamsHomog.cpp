#include "ConfigParamsHomog.h"
#include "ConfigFileReader.h"

bool ConfigParamsHomog::initParamsFromConfigFile(std::string &configFilePath)
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

		// homog parameters
		if ( !cfr.getTopLevelParameterValue("problem_specific", "input_file_path", homog.inputFilePath) )
		{
			return false;
		}
	}
	catch(...)
	{
		return false;
	}

	return true;
}