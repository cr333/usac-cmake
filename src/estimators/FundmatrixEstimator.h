#ifndef FUNDMATRIXESTIMATOR_H
#define FUNDMATRIXESTIMATOR_H

#include <iostream>
#include <fstream>
#include <string>
#include "config/ConfigParamsFundmatrix.h"
#include "utils/MathFunctions.h"
#include "utils/FundmatrixFunctions.h"
#include "utils/HomographyFunctions.h"
#include "USAC.h"

class FundMatrixEstimator: public USAC<FundMatrixEstimator>
{
	public:
		inline bool		 initProblem(const ConfigParamsFund& cfg, double* pointData);
		// ------------------------------------------------------------------------
		// storage for the final result
		std::vector<double> final_model_params_;
		std::vector<double> degen_final_model_params_;

	public:
		FundMatrixEstimator() 
		{
			input_points_ = NULL;
			data_matrix_  = NULL;
			degen_data_matrix_  = NULL;
			models_.clear();
			models_denorm_.clear();
		};
		~FundMatrixEstimator() 
		{
			if (input_points_) { delete[] input_points_; input_points_ = NULL; }
			if (data_matrix_) { delete[] data_matrix_; data_matrix_ = NULL; }
			if (degen_data_matrix_) { delete[] degen_data_matrix_; degen_data_matrix_ = NULL; }
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
		double*		 input_points_denorm_;					    // stores pointer to original input points
		double       degen_homog_threshold_;				    // threshold for h-degen test
		unsigned int degen_max_upgrade_samples_;				// maximum number of upgrade attempts
		USACConfig::MatrixDecomposition matrix_decomposition_method;  // QR/LU decomposition
		
		// ------------------------------------------------------------------------
		// temporary storage
		double* input_points_;								// stores normalized data points
		double* data_matrix_;								// linearized input data
		double* degen_data_matrix_;							// only for degeneracy testing
		std::vector<int> degen_outlier_flags_;				// for easy access to outliers to degeneracy
		double  m_T1_[9], m_T2_[9], m_T2_trans_[9];			// normalization matrices
		std::vector<double*> models_;						// stores vector of models
		std::vector<double*> models_denorm_;				// stores vector of (denormalized) models
};


// ============================================================================================
// initProblem: initializes problem specific data and parameters
// this function is called once per run on new data
// ============================================================================================
bool FundMatrixEstimator::initProblem(const ConfigParamsFund& cfg, double* pointData)
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
	MathTools::mattr(m_T2_trans_, m_T2_, 3, 3);

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
	data_matrix_ = new double[9*usac_num_data_points_];	// 9 values per correspondence
	FTools::computeDataMatrix(data_matrix_, usac_num_data_points_, input_points_);

	// if required, set up H-data matrix/storage for degeneracy testing
	degen_outlier_flags_.clear(); degen_outlier_flags_.resize(usac_num_data_points_, 0);
	if (usac_test_degeneracy_)
	{
		degen_homog_threshold_ = cfg.fund.hDegenThreshold;
		degen_max_upgrade_samples_ = cfg.fund.maxUpgradeSamples;
		degen_final_model_params_.clear(); degen_final_model_params_.resize(9);
		degen_data_matrix_ = new double[2*9*usac_num_data_points_];	// 2 equations per correspondence
		HTools::computeDataMatrix(degen_data_matrix_, usac_num_data_points_, input_points_);
	}
	else
	{
		degen_homog_threshold_ = 0.0;
		degen_max_upgrade_samples_ = 0;
		degen_final_model_params_.clear();
		degen_data_matrix_ = NULL;
	}

	// read in the f-matrix specific parameters from the config struct
	matrix_decomposition_method = cfg.fund.decompositionAlg;
	degen_homog_threshold_ = cfg.fund.hDegenThreshold*cfg.fund.hDegenThreshold;

	return true;
}


// ============================================================================================
// cleanupProblem: release any temporary problem specific data storage 
// this function is called at the end of each run on new data
// ============================================================================================
void FundMatrixEstimator::cleanupProblem()
{
	if (input_points_) { delete[] input_points_; input_points_ = NULL; }
	if (data_matrix_) { delete[] data_matrix_; data_matrix_ = NULL; }
	if (degen_data_matrix_) { delete[] degen_data_matrix_; degen_data_matrix_ = NULL; }
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
	degen_outlier_flags_.clear();
}


// ============================================================================================
// generateMinimalSampleModels: generates minimum sample model(s) from the data points whose  
// indices are currently stored in m_sample. 
// the generated models are stored in a vector of models and are all evaluated
// ============================================================================================
unsigned int FundMatrixEstimator::generateMinimalSampleModels()
{
	double A[9*9];
	unsigned int nsols = 0;

	// form the matrix of equations for this minimal sample
	double *src_ptr;
	double *dst_ptr = A;
	for (unsigned int i = 0; i < usac_min_sample_size_; ++i)
	{
		src_ptr = data_matrix_ + min_sample_[i];
		for (unsigned int j = 0; j < 9; ++j)
		{
			*dst_ptr = *src_ptr;
			++dst_ptr;
			src_ptr += usac_num_data_points_;
		}
	}

	// LU/QR factorization
	double sol[9*9];
	double poly[4], roots[3];
	double *f1, *f2;
	int nullbuff [18];
	f1 = sol;
	f2 = sol+9;
	if (matrix_decomposition_method == USACConfig::DECOMP_QR)
	{
		FTools::nullspaceQR7x9(A, sol);
	}
	else if (matrix_decomposition_method == USACConfig::DECOMP_LU)
	{
		for (unsigned int i = 7*9; i < 9*9; ++i)
		{
			A[i] = 0.0;
		}
		int nullsize = FTools::nullspace(A, f1, 9, nullbuff);
		if (nullsize != 2)
		{
			return 0;
		}
	}

	// solve polynomial
	FTools::makePolynomial(f1, f2, poly);  
	nsols = FTools::rroots3(poly, roots);

	// form up to three fundamental matrices
	double T2_F[9];
	for (unsigned int i = 0; i < nsols; ++i)
	{
		for (unsigned int j = 0; j < 9; ++j)
		{
			*(models_[i]+j) = f1[j] * roots[i] + f2[j] * (1 -roots[i]);
		}
		// store denormalized version as well
		MathTools::mmul(T2_F, m_T2_trans_, models_[i], 3);
		MathTools::mmul(models_denorm_[i], T2_F, m_T1_, 3);    
	}

	return nsols;
}


// ============================================================================================
// generateRefinedModel: compute model using non-minimal set of samples
// default operation is to use a weight of 1 for every data point
// ============================================================================================
bool FundMatrixEstimator::generateRefinedModel(std::vector<unsigned int>& sample,
											   const unsigned int numPoints,
											   bool weighted,
											   double* weights)
{
	// form the matrix of equations for this non-minimal sample
	double *A = new double[numPoints*9];	
	double *src_ptr;
	double *dst_ptr = A;
	for (unsigned int i = 0; i < numPoints; ++i)
	{
		src_ptr = data_matrix_ + sample[i];
		for (unsigned int j = 0; j < 9; ++j)
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
			src_ptr += usac_num_data_points_;
		}
	}

	double Cv[9*9];
	FTools::formCovMat(Cv, A, numPoints, 9);

	double V[9*9], D[9], *p;
	MathTools::svdu1v(D, Cv, 9, V, 9);

	unsigned int j = 0;
	for (unsigned int i = 1; i < 9; ++i)
	{
		if (D[i] < D[j]) 
		{
			j = i;
		}
	}
	p = V + j;

	for (unsigned int i = 0; i < 9; ++i)
	{
		*(models_[0]+i) = *p;
		p += 9;
	}
	FTools::singulF(models_[0]);
	// store denormalized version as well
	double T2_F[9];
	MathTools::mmul(T2_F, m_T2_trans_, models_[0], 3);
	MathTools::mmul(models_denorm_[0], T2_F, m_T1_, 3);   

	delete[] A;

	return true;
}


// ============================================================================================
// validateSample: check if minimal sample is valid
// here, just returns true
// ============================================================================================
bool FundMatrixEstimator::validateSample()
{
	return true;
}


// ============================================================================================
// validateModel: check if model computed from minimal sample is valid
// checks oriented constraints to determine model validity
// ============================================================================================
bool FundMatrixEstimator::validateModel(unsigned int modelIndex)
{
	// check oriented constraints
	double e[3], sig1, sig2;
	FTools::computeEpipole(e, models_[modelIndex]);

	sig1 = FTools::getOriSign(models_[modelIndex], e, input_points_ + 6*min_sample_[0]);
	for(unsigned int i = 1; i < min_sample_.size(); ++i)
	{
		sig2 = FTools::getOriSign(models_[modelIndex], e, input_points_ + 6*min_sample_[i]);
		if (sig1 * sig2 < 0) 
		{
			return false;
		}
	}
	return true;	
}


// ============================================================================================
// evaluateModel: test model against all/subset of the data points
// ============================================================================================
bool FundMatrixEstimator::evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested)
{
    double rx, ry, rwc, ryc, rxc, r, temp_err;
	double* model = models_denorm_[modelIndex];
	double* pt;
	std::vector<double>::iterator current_err_array = err_ptr_[0];
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

		// compute sampson error
		pt = input_points_denorm_ + 6*pt_index;
		rxc = (*model) * (*(pt+3)) + (*(model+3)) * (*(pt+4)) + (*(model+6));
		ryc = (*(model+1)) * (*(pt+3)) + (*(model+4)) * (*(pt+4)) + (*(model+7));
		rwc = (*(model+2)) * (*(pt+3)) + (*(model+5)) * (*(pt+4)) + (*(model+8));
		r =((*(pt)) * rxc + (*(pt+1)) * ryc + rwc);
		rx = (*model) * (*(pt)) + (*(model+1)) * (*(pt+1)) + (*(model+2));
		ry = (*(model+3)) * (*(pt)) + (*(model+4)) * (*(pt+1)) + (*(model+5)); 
		temp_err = r*r / (rxc*rxc + ryc*ryc + rx*rx + ry*ry);
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
// test if >=5 points in the sample are on a plane
// ============================================================================================
void FundMatrixEstimator::testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel)
{
	*degenerateModel = false;
	*upgradeModel = false;

	// make up the tuples to be used to check for degeneracy
	unsigned int degen_sample_indices[] = {0, 1, 2, 3,
										   3, 4, 5, 6,
										   0, 1, 5, 6,
										   0, 2, 4, 5,
										   1, 2, 4, 6,
										   0, 3, 4, 6,
										   1, 3, 4, 5,
										   2, 3, 5, 6};

	// the above tuples need to be tested on the remaining points for each case
	unsigned int test_point_indices[] = {4, 5, 6,
									     0, 1, 2,
									     2, 3, 4,
									     1, 3, 6,
									     0, 3, 5,
									     1, 2, 5,
									     0, 2, 6,
									     0, 1, 4};

	unsigned int *sample_pos = degen_sample_indices;
	unsigned int *test_pos = test_point_indices;
	double h[9];
	double T2_inv[9], T2_H[9];
	for (unsigned int i = 0; i < 9; ++i)
	{
		T2_inv[i] = m_T2_[i];
	}
	MathTools::minv(T2_inv, 3);

	std::vector<unsigned int> sample(7), test(3);
	std::vector<double> errs;
	for(unsigned int i = 0; i < 8; ++i)
	{
		// compute H from the current set of 4 points
		for (unsigned int j = 0; j < 4; ++j)
		{
			sample[j] = min_sample_[sample_pos[j]];
		}
		FTools::computeHFromMinCorrs(sample, 4, usac_num_data_points_, degen_data_matrix_, h);
		MathTools::mmul(T2_H, T2_inv, h, 3);
		MathTools::mmul(h, T2_H, m_T1_, 3);  

		// check test points to see how many are consistent
		for (unsigned int j = 0; j < 3; ++j)
		{
			test[j] = min_sample_[test_pos[j]];
		}
		unsigned int num_inliers = FTools::getHError(test, 3, errs, input_points_denorm_, h, degen_homog_threshold_);
		for (unsigned int j = 0, count = 4; j < 3; ++j)
		{
			if (errs[j] < degen_homog_threshold_)
			{
				sample[count++] = test[j];
			}
		}

		// if at least 1 inlier in the test points, then h-degenerate sample found
		if (num_inliers > 0)
		{
			// find inliers from all data points
			num_inliers = FTools::getHError(evaluation_pool_, usac_num_data_points_, errs, input_points_denorm_, h, degen_homog_threshold_);
			//std::cout << "Degenerate sample found with " << num_inliers << " inliers" << std::endl;

			// refine with least squares fit
			unsigned int count = 0;
			std::vector<unsigned int> inlier_sample(num_inliers);
			for (unsigned int j = 0; j < usac_num_data_points_; ++j)
			{
				if (errs[j] < degen_homog_threshold_)
				{
					inlier_sample[count++] = evaluation_pool_[j];
				}
			}
			FTools::computeHFromCorrs(inlier_sample, inlier_sample.size(), usac_num_data_points_, degen_data_matrix_, h);
			MathTools::mmul(T2_H, T2_inv, h, 3);
			MathTools::mmul(h, T2_H, m_T1_, 3);  

			// find support of homography
			num_inliers = FTools::getHError(evaluation_pool_, usac_num_data_points_, errs, input_points_denorm_, h, degen_homog_threshold_);
			//std::cout << "Refined model has " << num_inliers << " inliers" << std::endl;
#if 1
			if (num_inliers < usac_results_.best_inlier_count_/5)
			{
				sample_pos += 4;
				test_pos += 3;
				continue;
			}
#endif
			// set flag
			*degenerateModel = true;

			// if largest degenerate model found so far, store results
			if (num_inliers > usac_results_.degen_inlier_count_)
			{
				// set flag
				*upgradeModel = true;

				// refine with final least squares fit
				count = 0;
				inlier_sample.resize(num_inliers);
				for (unsigned int j = 0; j < usac_num_data_points_; ++j)
				{
					if (errs[j] < degen_homog_threshold_)
					{
						inlier_sample[count++] = evaluation_pool_[j];
					}
				}
				FTools::computeHFromCorrs(inlier_sample, inlier_sample.size(), usac_num_data_points_, degen_data_matrix_, h);

				usac_results_.degen_inlier_count_ = num_inliers;
				// store homography
				for (unsigned int j = 0; j < 9; ++j)
				{
					degen_final_model_params_[j] = h[j];
				}
				// store inliers and outliers - for use in model completion
				for (unsigned int j = 0; j < usac_num_data_points_; ++j)
				{
					if (errs[j] < degen_homog_threshold_)
					{
						usac_results_.degen_inlier_flags_[evaluation_pool_[j]] = 1;
						degen_outlier_flags_[evaluation_pool_[j]] = 0;
					}
					else
					{
						degen_outlier_flags_[evaluation_pool_[j]] = 1;
						usac_results_.degen_inlier_flags_[evaluation_pool_[j]] = 0;
					}
				}
				// store the degenerate points from the minimal sample
				usac_results_.degen_sample_ = sample;

			} // end store denerate results
		} // end check for one model degeneracy

		sample_pos += 4;
		test_pos += 3;

	} // end check for all combinations in the minimal sample
}


// ============================================================================================
// upgradeDegenerateModel: try to upgrade degenerate model to non-degenerate by sampling from
// the set of outliers to the degenerate model
// ============================================================================================
unsigned int FundMatrixEstimator::upgradeDegenerateModel()
{
	unsigned int best_upgrade_inliers = usac_results_.best_inlier_count_;
	unsigned int num_outliers = usac_num_data_points_ - usac_results_.degen_inlier_count_;

	if (num_outliers < 2) {
		return 0;
	}

	std::vector<unsigned int> outlier_indices(num_outliers);
	unsigned int count = 0;
	for (unsigned int i = 0; i < usac_num_data_points_; ++i)
	{
		if (degen_outlier_flags_[i])
		{
			outlier_indices[count++] = i;
		}
	}
	std::vector<unsigned int> outlier_sample(2);
	std::vector<double>::iterator current_err_array = err_ptr_[0];

	double* pt1_index, *pt2_index;
	double x1[3], x1p[3], x2[3], x2p[3];
	double temp[3], l1[3], l2[3], ep[3];
	double skew_sym_ep[9];
	double T2_F[9];
	for (unsigned int i = 0; i < degen_max_upgrade_samples_; ++i)
	{
		generateUniformRandomSample(num_outliers, 2, &outlier_sample);

		pt1_index = input_points_ + 6*outlier_indices[outlier_sample[0]];
		pt2_index = input_points_ + 6*outlier_indices[outlier_sample[1]];

		x1[0]  = pt1_index[0]; x1[1]  = pt1_index[1]; x1[2]  = 1.0;
		x1p[0] = pt1_index[3]; x1p[1] = pt1_index[4]; x1p[2] = 1.0;
		x2[0]  = pt2_index[0]; x2[1]  = pt2_index[1]; x2[2]  = 1.0;
		x2p[0] = pt2_index[3]; x2p[1] = pt2_index[4]; x2p[2] = 1.0;

		MathTools::vmul(temp, &degen_final_model_params_[0], x1, 3);
		MathTools::crossprod(l1, temp, x1p, 1);

		MathTools::vmul(temp, &degen_final_model_params_[0], x2, 3);
		MathTools::crossprod(l2, temp, x2p, 1);

		MathTools::crossprod(ep, l1, l2, 1);

		MathTools::skew_sym(skew_sym_ep, ep);
		MathTools::mmul(models_[0], skew_sym_ep, &degen_final_model_params_[0], 3);
		MathTools::mmul(T2_F, m_T2_trans_, models_[0], 3);
		MathTools::mmul(models_denorm_[0], T2_F, m_T1_, 3);   

		unsigned int num_inliers, num_pts_tested;
		evaluateModel(0, &num_inliers, &num_pts_tested);

		if (num_inliers > best_upgrade_inliers)
		{
			usac_results_.degen_sample_[5] = outlier_indices[outlier_sample[0]];
			usac_results_.degen_sample_[6] = outlier_indices[outlier_sample[1]];
			min_sample_ = usac_results_.degen_sample_;
			storeSolution(0, num_inliers);
			best_upgrade_inliers = num_inliers;

			unsigned int count = 0;
			for (size_t j = 0; j < outlier_indices.size(); ++j)
			{
				if (*(current_err_array+outlier_indices[j]) < usac_inlier_threshold_)
				{
					++count;
				}
			}
			unsigned int num_samples = updateStandardStopping(count, num_outliers, 2);
			//std::cout << "Inliers = " << num_inliers << ", in/out = " << count << "/" << num_outliers 
			//	      << ". Num samples = " << num_samples << std::endl;
			if (num_samples < degen_max_upgrade_samples_)
			{
				degen_max_upgrade_samples_ = num_samples;
			}

		}
	}

	std::cout << "Upgraded model has " << best_upgrade_inliers << " inliers" << std::endl;
	return best_upgrade_inliers;
}


// ============================================================================================
// findWeights: given model and points, compute weights to be used in local optimization
// ============================================================================================
void FundMatrixEstimator::findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers, 
									  unsigned int numInliers, double* weights)
{
    double rx, ry, ryc, rxc;
	double* model = models_[modelIndex];
	double* pt;
	unsigned int pt_index;

	for (unsigned int i = 0; i < numInliers; ++i)
	{
		// get index of point to be verified
		pt_index = inliers[i];

		// compute weight (ref: torr dissertation, eqn. 2.25)
		pt = input_points_ + 6*pt_index;
		rxc = (*model) * (*(pt+3)) + (*(model+3)) * (*(pt+4)) + (*(model+6));
		ryc = (*(model+1)) * (*(pt+3)) + (*(model+4)) * (*(pt+4)) + (*(model+7));
		rx = (*model) * (*(pt)) + (*(model+1)) * (*(pt+1)) + (*(model+2));
		ry = (*(model+3)) * (*(pt)) + (*(model+4)) * (*(pt+1)) + (*(model+5)); 

		weights[i] = 1/sqrt(rxc*rxc + ryc*ryc + rx*rx + ry*ry);
	}
}


// ============================================================================================
// storeModel: stores current best model
// this function is called  (by USAC) every time a new best model is found
// ============================================================================================
void FundMatrixEstimator::storeModel(unsigned int modelIndex, unsigned int numInliers)
{
	// save the current model as the best solution so far
	for (unsigned int i = 0; i < 9; ++i)
	{
		final_model_params_[i] = *(models_denorm_[modelIndex]+i);
	}
}

#endif

