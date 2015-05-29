#ifndef TEMPLATEESTIMATOR_H
#define TEMPLATEESTIMATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include "config/ConfigParamsTemplate.h"
#include "utils/MathFunctions.h"
#include "USAC.h"

class TemplateEstimator: public USAC<TemplateEstimator>
{
	public:
		inline bool		 initProblem(const ConfigParamsTemplate& cfg, double* pointData);

		// ------------------------------------------------------------------------
		// problem specific functions: implement these
		inline void		 cleanupProblem();
		inline unsigned int generateMinimalSampleModels();
		inline bool		 generateRefinedModel(std::vector<unsigned int>& sample, const unsigned int numPoints, 
										  bool weighted = false, double* weights = NULL);
		inline bool		 validateSample();
		inline bool		 validateModel(unsigned int modelIndex);
		inline bool		 evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested);
		inline void		 testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel);
		inline unsigned int upgradeDegenerateModel();
		inline void		 findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers, 
									 unsigned int numInliers, double* weights);
		inline void		 storeModel(unsigned int modelIndex, unsigned int numInliers);

	public:
		// ------------------------------------------------------------------------
		// storage for the final result
		std::vector<double> final_model_params_;

	private:
};


// ============================================================================================
// initProblem: initializes problem specific data and parameters
// call this function once per run on new data
// ============================================================================================
bool TemplateEstimator::initProblem(const ConfigParamsTemplate& cfg, double* pointData)
{
	// set up problem specific data here (e.g., data matrices, parameter vectors, etc.)
	return true;
}


// ============================================================================================
// cleanupProblem: release any temporary problem specific data storage 
// call this function once at the end of each run on new data
// ============================================================================================
void TemplateEstimator::cleanupProblem()
{

}


// ============================================================================================
// generateMinimalSampleModels: generates minimum sample model from the data points whose  
// indices are currently stored in min_sample_. 
// in this case, only one model per minimum sample
// ============================================================================================
unsigned int TemplateEstimator::generateMinimalSampleModels()
{
	// this function must be implemented
	
	return num_models;  // return the number of minimal sample models
}


// ============================================================================================
// generateRefinedModel: compute model using non-minimal set of samples
// default operation is to use a weight of 1 for every data point
// ============================================================================================
bool TemplateEstimator::generateRefinedModel(std::vector<unsigned int>& sample,
										  unsigned int numPoints,
										  bool weighted,
										  double* weights)
{
	// implement this function if you want to use local optimization

	return true;
}


// ============================================================================================
// validateSample: check if minimal sample is valid
// ============================================================================================
bool TemplateEstimator::validateSample()
{
   // implement this function if you want to pre-verify the minimal sample, otherwise
   // simply return true
   return true;	
}


// ============================================================================================
// validateModel: check if model computed from minimal sample is valid
// ============================================================================================
bool TemplateEstimator::validateModel(const unsigned int modelIndex)
{
   // implement this function if you want to pre-verify a model, otherwise
   // simply return true
	return true;
}


// ============================================================================================
// evaluateModel: test model against all/subset of the data points
// ============================================================================================
bool TemplateEstimator::evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested)
{
	// test model against all data points, or a randomized subset (in the case of early
	// termination with, for e.g., SPRT)
	double* model = models_denorm_[modelIndex];
	std::vector<double>::iterator current_err_array = err_ptr_[0];
	double err;
	bool good_flag = true;
	double lambdaj, lambdaj_1 = 1.0;
	*numInliers = 0;
	*numPointsTested = 0;
	unsigned int pt_index;

	for (unsigned int i = 0; i < usac_num_data_points_; ++i)
	{
		// get index of point to be verified
		if (eval_pool_index_ > usac_num_data_points_-1)
		{
			eval_pool_index_ = 0;
		}
		pt_index = evaluation_pool_[eval_pool_index_];
		++eval_pool_index_;

		// compute point-model error for the problem of interest
		//
		// --- implement this
		//

		*(current_err_array+pt_index) = err;

		if (err < usac_inlier_threshold_)
		{
			++(*numInliers);
		}

		if (usac_verif_method_ == USACConfig::VERIF_SPRT)
		{
			if (err < usac_inlier_threshold_)
			{			
				lambdaj = lambdaj_1 * (sprt_delta_/sprt_epsilon_);
			}
			else
			{
				lambdaj = lambdaj_1 * ( (1 - sprt_delta_)/(1 - sprt_epsilon_) );
			}

			if (lambdaj > decision_threshold_sprt_)
			{
				good_flag = false;
				*numPointsTested = i+1;
				return good_flag;
			}
			else
			{
				lambdaj_1 = lambdaj;
			}
		}
	}
	*numPointsTested = usac_num_data_points_;

	return good_flag;
}

// ============================================================================================
// testSolutionDegeneracy: check if model is degenerate
// ============================================================================================
void TemplateEstimator::testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel)
{
	// implement this if you want to detect degenerate models, otherwise return false
	*degenerateModel = false;
	*upgradeModel = false;
}

// ============================================================================================
// upgradeDegenerateModel: try to upgrade degenerate model to non-degenerate by sampling from
// the set of outliers to the degenerate model
// ============================================================================================
unsigned int TemplateEstimator::upgradeDegenerateModel()
{
	// implement this if you want to upgrade degenerate models to non-degenerate, otherwise
	// return 0;
	return num_inliers;
}


// ============================================================================================
// findWeights: given model and points, compute weights to be used in local optimization
// ============================================================================================
void TemplateEstimator::findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers, 
								 unsigned int numInliers, double* weights)
{
	// implement this if you want weighted local refinement, otherwise return an array of ones
	for (unsigned int i = 0; i < numInliers; ++i)
	{
		weights[i] = 1.0;
	}
}


// ============================================================================================
// storeModel: stores current best model
// this function is called  (by USAC) every time a new best model is found
// ============================================================================================
void TemplateEstimator::storeModel(const unsigned int modelIndex, unsigned int numInliers)
{
	// save the current model as the best solution so far
}

#endif

