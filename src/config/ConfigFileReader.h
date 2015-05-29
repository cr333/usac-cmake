#ifndef CONFIGFILEREADER_H
#define CONFIGFILEREADER_H
#pragma warning(disable: 4290)	// disable libconfig compile warnings
#include <cstdlib>
#include <string>
#include <iostream>
#include <libconfig.h++>

class ConfigFileReader: public libconfig::Config
{
    public:
		template <class T>
		bool getTopLevelParameterValue(const std::string settingName, const std::string parameterName, 
									   T& value) const
		{
			const libconfig::Setting& root = getRoot();
			try
			{
				const libconfig::Setting &setting = root[settingName];
				if (!setting.lookupValue(parameterName.c_str(), value))
				{
					std::cout << settingName << ": parameter " << parameterName << " not found" << std::endl;
					return false;
				}
			} // end try 
			catch(...)
			{
				std::cerr << settingName << " block not found; recheck configuration file" << std::endl;
				return false;
			}
			return true;
		}

		template <class T>
		bool getSubLevelParameterValue(const std::string settingName, const std::string subSettingName,
									   const std::string parameterName, T& value) const
		{
			const libconfig::Setting& root = getRoot();
			try
			{
				const libconfig::Setting &setting = root[settingName][subSettingName];
				if (!setting.lookupValue(parameterName.c_str(), value))
				{
					std::cout << settingName << ":" << subSettingName
							  << ": parameter " << parameterName << " not found" << std::endl;
					return false;
				}
			} 
			catch(...)
			{
				std::cerr << settingName << ":" << subSettingName << " block not found; recheck configuration file" << std::endl;
				return false;
			}
			return true;
		}

};

#endif