#include "ConfigFileReader.h"
#include "ConfigParams.h"

bool ConfigParams::initUSACParamsFromConfigFile(const ConfigFileReader& cfr)
{
	// read in parameters 
	try
	{
		// first get all common parameters
		std::string rand_sampling_method, verif_method, local_opt_method, est_problem;
		if( !( cfr.getTopLevelParameterValue("common", "ransac_conf", common.confThreshold) &&
			   cfr.getTopLevelParameterValue("common", "min_sample_size", common.minSampleSize) &&
			   cfr.getTopLevelParameterValue("common", "inlier_threshold", common.inlierThreshold) &&
			   cfr.getTopLevelParameterValue("common", "max_hypotheses", common.maxHypotheses) &&
			   cfr.getTopLevelParameterValue("common", "max_solutions_per_sample", common.maxSolutionsPerSample) &&
			   cfr.getTopLevelParameterValue("common", "prevalidate_sample", common.prevalidateSample) &&
			   cfr.getTopLevelParameterValue("common", "prevalidate_model", common.prevalidateModel) &&
			   cfr.getTopLevelParameterValue("common", "degen_check", common.testDegeneracy) &&
			   cfr.getTopLevelParameterValue("common", "random_sampling_method", rand_sampling_method) &&
			   cfr.getTopLevelParameterValue("common", "verif_method", verif_method) &&
			   cfr.getTopLevelParameterValue("common", "local_opt", local_opt_method) ) )
		{
			return false;
		}
		else
		{
			// verify parameter values 
			if (common.confThreshold < 0 || common.confThreshold > 1)
			{
				std::cout << "RANSAC confidence value must be between 0 and 1" << std::endl;
				return false;
			}

			if (rand_sampling_method.compare("UNIFORM") == 0)
				common.randomSamplingMethod = USACConfig::SAMP_UNIFORM;
			else if (rand_sampling_method.compare("PROSAC") == 0)
				common.randomSamplingMethod = USACConfig::SAMP_PROSAC;
			else
			{
				std::cerr << "Random sampling method " << rand_sampling_method << " not recognized" << std::endl;
				return false;
			}

			if (verif_method.compare("STANDARD") == 0)
				common.verifMethod = USACConfig::VERIF_STANDARD;
			else if (verif_method.compare("SPRT") == 0)
				common.verifMethod = USACConfig::VERIF_SPRT;
			else
			{
				std::cerr << "Verification method " << verif_method << " not recognized" << std::endl;
				return false;
			}

			if (local_opt_method.compare("NO_LO") == 0)
				common.localOptMethod = USACConfig::LO_NONE;
			else if (local_opt_method.compare("LOSAC") == 0)
				common.localOptMethod = USACConfig::LO_LOSAC;
			else
			{
				std::cerr << "Local optimization method " << local_opt_method << " not recognized" << std::endl;
				return false;
			}
		} 

		// read in PROSAC parameters if required
		if (common.randomSamplingMethod == USACConfig::SAMP_PROSAC)
		{
			if (!( cfr.getTopLevelParameterValue("prosac", "max_prosac_samples", prosac.maxSamples) &&
				   cfr.getTopLevelParameterValue("prosac", "beta", prosac.beta) &&
				   cfr.getTopLevelParameterValue("prosac", "non_rand_conf", prosac.nonRandConf) &&
				   cfr.getTopLevelParameterValue("prosac", "min_stopping_length", prosac.minStopLen) &&
			       cfr.getTopLevelParameterValue("prosac", "sorted_points_path", prosac.sortedPointsFile) ))
			{
				return false;
			}
		}

		// read in SPRT parameters if required
		if (common.verifMethod == USACConfig::VERIF_SPRT)
		{
			if ( !(cfr.getTopLevelParameterValue("sprt", "time_model", sprt.tM) &&
			       cfr.getTopLevelParameterValue("sprt", "models_per_sample", sprt.mS) &&
			       cfr.getTopLevelParameterValue("sprt", "delta", sprt.delta) &&
			       cfr.getTopLevelParameterValue("sprt", "epsilon", sprt.epsilon)) )
			{
				return false;
			}
		}

		// read in local optimization parameters if required
		if (common.localOptMethod == USACConfig::LO_LOSAC)
		{
			if ( !(cfr.getTopLevelParameterValue("losac", "inner_sample_size", losac.innerSampleSize) &&
				   cfr.getTopLevelParameterValue("losac", "inner_ransac_repetitions", losac.innerRansacRepetitions) &&
				   cfr.getTopLevelParameterValue("losac", "threshold_multiplier", losac.thresholdMultiplier) &&
				   cfr.getTopLevelParameterValue("losac", "num_steps", losac.numStepsIterative)) )
			{
				return false;
			}
		}	
	}
	catch(...)
	{
		return false;
	}

	return true;
}

bool ConfigParams::initParamsFromConfigFile(std::string& configFilePath) {
	std::cerr << "Implement this in the derived class" << std::endl;
	return false;
}