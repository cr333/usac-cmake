#ifndef HOMOGESTIMATOR_H
#define HOMOGESTIMATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include "config/ConfigParamsHomog.h"
#include "utils/MathFunctions.h"
#include "utils/FundmatrixFunctions.h"
#include "utils/HomographyFunctions.h"
#include "USAC.h"

class HomogEstimator: public USAC<HomogEstimator>
{
	public:
		inline bool		 initProblem(const ConfigParamsHomog& cfg, double* pointData);
		// ------------------------------------------------------------------------
		// storage for the final result
		std::vector<double> final_model_params_;

	public:
		HomogEstimator() 
		{
			input_points_ = NULL;
			data_matrix_  = NULL;
			models_.clear();
			models_denorm_.clear();
		};
		~HomogEstimator() 
		{
			if (input_points_) { delete[] input_points_; input_points_ = NULL; }
			if (data_matrix_) { delete[] data_matrix_; data_matrix_ = NULL; }
			for (size_t i = 0; i < models_.size(); ++i)
			{
				if (models_[i]) { delete[] models_[i]; }
			}
			models_.clear();
			for (size_t i = 0; i < models_denorm_.size(); ++i)
			{
				if (models_denorm_[i]) { delete[] models_denorm_[i]; }
			}
			models_denorm_.clear();
		};

	public:
		// ------------------------------------------------------------------------
		// problem specific functions
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

	private:
		double*      input_points_denorm_;					// stores pointer to original input points

		// ------------------------------------------------------------------------
		// temporary storage
		double* input_points_;							// stores normalized data points
		double* data_matrix_;							// linearized input data
		double  m_T1_[9], m_T2_[9], m_T2inv_[9];			// normalization matrices
		std::vector<double*> models_;				    // stores vector of models
		std::vector<double*> models_denorm_;			// stores vector of (denormalized) models
};


// ============================================================================================
// initProblem: initializes problem specific data and parameters
// this function is called once per run on new data
// ============================================================================================
bool HomogEstimator::initProblem(const ConfigParamsHomog& cfg, double* pointData)
{
	// copy pointer to input data
	input_points_denorm_ = pointData;
	input_points_       = new double[6*cfg.common.numDataPoints];
	if (input_points_denorm_ == NULL)
	{
		std::cerr << "Input point data not properly initialized" << std::endl;
		return false;
	}
	if (input_points_ == NULL)
	{
		std::cerr << "Could not allocate storage for normalized data points" << std::endl;
		return false;
	}

	// normalize input data
	// following this, input_points_ has the normalized points and input_points_denorm_ has 
	// the original input points
	FTools::normalizePoints(input_points_denorm_, input_points_, cfg.common.numDataPoints, m_T1_, m_T2_);
	for (unsigned int i = 0; i < 9; ++i)
	{
		m_T2inv_[i] = m_T2_[i];
	}
	MathTools::minv(m_T2inv_, 3);

	// allocate storage for models
	final_model_params_.clear(); final_model_params_.resize(9);
	models_.clear(); models_.resize(usac_max_solns_per_sample_);
	models_denorm_.clear(); models_denorm_.resize(usac_max_solns_per_sample_);
	for (unsigned int i = 0; i < usac_max_solns_per_sample_; ++i)
	{
		models_[i] = new double[9];
		models_denorm_[i] = new double[9];
	}

	// precompute the data matrix
	data_matrix_ = new double[18*usac_num_data_points_];	// 2 equations per correspondence
	HTools::computeDataMatrix(data_matrix_, usac_num_data_points_, input_points_);

	return true;
}


// ============================================================================================
// cleanupProblem: release any temporary problem specific data storage 
// this function is called at the end of each run on new data
// ============================================================================================
void HomogEstimator::cleanupProblem()
{
	if (input_points_) { delete[] input_points_; input_points_ = NULL; }
	if (data_matrix_) { delete[] data_matrix_; data_matrix_ = NULL; }
	for (size_t i = 0; i < models_.size(); ++i)
	{
		if (models_[i]) { delete[] models_[i]; }
	}
	models_.clear();
	for (size_t i = 0; i < models_denorm_.size(); ++i)
	{
		if (models_denorm_[i]) { delete[] models_denorm_[i]; }
	}
	models_denorm_.clear();
}


// ============================================================================================
// generateMinimalSampleModels: generates minimum sample model from the data points whose  
// indices are currently stored in min_sample_. 
// in this case, only one model per minimum sample
// ============================================================================================
unsigned int HomogEstimator::generateMinimalSampleModels()
{
   double A[8*9];
   double At[9*8];

	// form the matrix of equations for this minimal sample
	double *src_ptr;
	double *dst_ptr = A;
	for (unsigned int i = 0; i < usac_min_sample_size_; ++i)
	{
		for (unsigned int j = 0; j < 2; ++j)
		{
			src_ptr = data_matrix_ + 2*min_sample_[i] + j;
			for (unsigned int k = 0; k < 9; ++k)
			{
				*dst_ptr = *src_ptr; 
				++dst_ptr;
				src_ptr += 2*usac_num_data_points_;
			}
		}
	}

	MathTools::mattr(At, A, 8, 9);

	double D[9], U[9*9], V[8*8], *p;
	MathTools::svduv(D, At, U, 9, V, 8);
	p = U + 8;

	double T2_H[9];
	for (unsigned int i = 0; i < 9; ++i)
	{
		*(models_[0]+i) = *p;
		p += 9;
	}
	MathTools::mmul(T2_H, m_T2inv_, models_[0], 3);
	MathTools::mmul(models_denorm_[0], T2_H, m_T1_, 3);  

	return 1;
}


// ============================================================================================
// generateRefinedModel: compute model using non-minimal set of samples
// default operation is to use a weight of 1 for every data point
// ============================================================================================
bool HomogEstimator::generateRefinedModel(std::vector<unsigned int>& sample,
										  unsigned int numPoints,
										  bool weighted,
										  double* weights)
{
	// form the matrix of equations for this non-minimal sample
	double *A = new double[numPoints*2*9];	
	double *src_ptr;
	double *dst_ptr = A;
	for (unsigned int i = 0; i < numPoints; ++i)
	{
		for (unsigned int j = 0; j < 2; ++j)
		{
			src_ptr = data_matrix_ + 2*sample[i] + j;
			for (unsigned int k = 0; k < 9; ++k)
			{
				if (!weighted)
				{
					*dst_ptr = *src_ptr;
				}
				else
				{
					*dst_ptr = (*src_ptr)*weights[i];
				}
				++dst_ptr;
				src_ptr += 2*usac_num_data_points_;
			}
		}
	}

	// decompose
	double V[9*9], D[9], *p;
	MathTools::svdu1v(D, A, 2*numPoints, V, 9);

	unsigned int j = 0;
	for (unsigned int i = 1; i < 9; ++i)
	{
		if (D[i] < D[j]) 
			j = i;
	}
	p = V + j;

	for (unsigned int i = 0; i < 9; ++i)
	{
		*(models_[0]+i) = *p;
		p += 9;
	}
	double T2_H[9];
	MathTools::mmul(T2_H, m_T2inv_, models_[0], 3);
	MathTools::mmul(models_denorm_[0], T2_H, m_T1_, 3); 

	delete[] A;

	return true;
}


// ============================================================================================
// validateSample: check if minimal sample is valid
// ============================================================================================
bool HomogEstimator::validateSample()
{
	// check oriented constraints
   double p[3], q[3];
   double *a, *b, *c, *d;

   a = input_points_ + 6*min_sample_[0];
   b = input_points_ + 6*min_sample_[1];
   c = input_points_ + 6*min_sample_[2];
   d = input_points_ + 6*min_sample_[3];

   HTools::crossprod(p, a, b, 1);
   HTools::crossprod(q, a+3, b+3, 1);

   if ((p[0]*c[0]+p[1]*c[1]+p[2]*c[2])*(q[0]*c[3]+q[1]*c[4]+q[2]*c[5])<0)
      return false;
   if ((p[0]*d[0]+p[1]*d[1]+p[2]*d[2])*(q[0]*d[3]+q[1]*d[4]+q[2]*d[5])<0)
      return false;

   HTools::crossprod(p, c, d, 1);
   HTools::crossprod(q, c+3, d+3, 1);

   if ((p[0]*a[0]+p[1]*a[1]+p[2]*a[2])*(q[0]*a[3]+q[1]*a[4]+q[2]*a[5])<0)
      return false;
   if ((p[0]*b[0]+p[1]*b[1]+p[2]*b[2])*(q[0]*b[3]+q[1]*b[4]+q[2]*b[5])<0)
      return false;

   return true;	
}


// ============================================================================================
// validateModel: check if model computed from minimal sample is valid
// ============================================================================================
bool HomogEstimator::validateModel(const unsigned int modelIndex)
{
	return true;
}


// ============================================================================================
// evaluateModel: test model against all/subset of the data points
// ============================================================================================
bool HomogEstimator::evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested)
{
	double* model = models_denorm_[modelIndex];
	double inv_model[9];
	double h_x[3], h_inv_xp[3], temp_err;
	double* pt;
	std::vector<double>::iterator current_err_array = err_ptr_[0];
	bool good_flag = true;
	double lambdaj, lambdaj_1 = 1.0;
	*numInliers = 0;
	*numPointsTested = 0;
	unsigned int pt_index;

	for (unsigned int i = 0; i < 9; ++i)
	{
		inv_model[i] = model[i];
	}
	MathTools::minv(inv_model, 3);
	for (unsigned int i = 0; i < usac_num_data_points_; ++i)
	{
		// get index of point to be verified
		if (eval_pool_index_ > usac_num_data_points_-1)
		{
			eval_pool_index_ = 0;
		}
		pt_index = evaluation_pool_[eval_pool_index_];
		++eval_pool_index_;

		// compute symmetric transfer error
		pt = input_points_denorm_ + 6*pt_index;
		MathTools::vmul(h_x, model, pt, 3);
		MathTools::vmul(h_inv_xp, inv_model, pt+3, 3);

		double err1 = 0.0, err2 = 0.0;
		for (unsigned int j = 0; j < 2; ++j)
		{
			err1 += (h_x[j]/h_x[2] - pt[3+j]) * (h_x[j]/h_x[2] - pt[3+j]);
			err2 += (h_inv_xp[j]/h_inv_xp[2] - pt[j]) * (h_inv_xp[j]/h_inv_xp[2] - pt[j]);
		}
		temp_err = err1 + err2;
		*(current_err_array+pt_index) = temp_err;

		if (temp_err < usac_inlier_threshold_)
		{
			++(*numInliers);
		}

		if (usac_verif_method_ == USACConfig::VERIF_SPRT)
		{
			if (temp_err < usac_inlier_threshold_)
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
void HomogEstimator::testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel)
{
	*degenerateModel = false;
	*upgradeModel = false;
}

// ============================================================================================
// upgradeDegenerateModel: try to upgrade degenerate model to non-degenerate by sampling from
// the set of outliers to the degenerate model
// ============================================================================================
unsigned int HomogEstimator::upgradeDegenerateModel()
{
	return 0;
}


// ============================================================================================
// findWeights: given model and points, compute weights to be used in local optimization
// ============================================================================================
void HomogEstimator::findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers, 
								 unsigned int numInliers, double* weights)
{
	for (unsigned int i = 0; i < numInliers; ++i)
	{
		weights[i] = 1.0;
	}
}


// ============================================================================================
// storeModel: stores current best model
// this function is called  (by USAC) every time a new best model is found
// ============================================================================================
void HomogEstimator::storeModel(const unsigned int modelIndex, unsigned int numInliers)
{
	// save the current model as the best solution so far
	for (unsigned int i = 0; i < 9; ++i)
	{
		final_model_params_[i] = *(models_denorm_[modelIndex]+i);
	}
}

#endif

