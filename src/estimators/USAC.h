#ifndef USAC_HH
#define USAC_HH

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#endif

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>
#include "config/ConfigParams.h"
#include "utils/Timer.h"

struct UsacResults {
	void reset() {
		hyp_count_			     = 0;
		model_count_		     = 0;
		rejected_sample_count_   = 0;
		rejected_model_count_    = 0;
		best_inlier_count_       = 0;
		degen_inlier_count_		 = 0;
		total_points_verified_   = 0;
		num_local_optimizations_ = 0;
		total_runtime_			 = 0;
		inlier_flags_.clear();
		best_sample_.clear();
		degen_inlier_flags_.clear();
		degen_sample_.clear();
	}
	unsigned int			  hyp_count_;
	unsigned int			  model_count_;
	unsigned int			  rejected_sample_count_;
	unsigned int			  rejected_model_count_;
	unsigned int			  total_points_verified_;
	unsigned int			  num_local_optimizations_;
	unsigned int			  best_inlier_count_;
	std::vector<unsigned int> inlier_flags_;
	std::vector<unsigned int> best_sample_;
	unsigned int			  degen_inlier_count_;
	std::vector<unsigned int> degen_inlier_flags_;
	std::vector<unsigned int> degen_sample_;
	double					  total_runtime_;
};

template <class ProblemType>
class USAC
{
	public:
		// ------------------------------------------------------------------------
		//  initialization/main functions
		USAC() {};
		virtual ~USAC() {};
		void initParamsUSAC(const ConfigParams& cfg);
		void initDataUSAC(const ConfigParams& cfg);
		bool solve();
		UsacResults usac_results_;
		
	protected:
		// ------------------------------------------------------------------------
		// pure virtual functions to be overridden in the derived class
		virtual inline unsigned int generateMinimalSampleModels() = 0;
		virtual inline bool		    generateRefinedModel(std::vector<unsigned int>& sample, const unsigned int numPoints, 
														 bool weighted, double* weights) = 0;
		virtual inline bool		    validateSample() = 0;
		virtual inline bool		    validateModel(unsigned int modelIndex) = 0;
		virtual inline bool		    evaluateModel(unsigned int modelIndex, unsigned int* numInliers, unsigned int* numPointsTested) = 0;
		virtual inline void		    testSolutionDegeneracy(bool* degenerateModel, bool* upgradeModel) = 0;
		virtual inline unsigned int upgradeDegenerateModel() = 0;
		virtual inline void		    findWeights(unsigned int modelIndex, const std::vector<unsigned int>& inliers, 
											    unsigned int numInliers, double* weights) = 0;
		virtual inline void		    storeModel(unsigned int modelIndex, unsigned int numInliers) = 0;

	// ------------------------------------------------------------------------
	// USAC input parameters
	protected:
		// common parameters
		double			usac_conf_threshold_;
		unsigned int    usac_min_sample_size_;  
		double			usac_inlier_threshold_;   
		unsigned int	usac_max_hypotheses_;	
		unsigned int	usac_max_solns_per_sample_; 
		bool			usac_prevalidate_sample_;
		bool			usac_prevalidate_model_;
		bool			usac_test_degeneracy_;
		unsigned int	usac_num_data_points_;  
		USACConfig::RandomSamplingMethod	   usac_sampling_method_;
		USACConfig::VerifMethod				   usac_verif_method_;
		USACConfig::LocalOptimizationMethod    usac_local_opt_method_;

		// PROSAC parameters
		unsigned int  prosac_growth_max_samples_;
		double        prosac_beta_;
		double		  prosac_non_rand_conf_;
		unsigned int  prosac_min_stop_length_;
		std::vector<int>   prosac_sorted_point_indices_;

		// SPRT parameters
		double sprt_tM_;
		double sprt_mS_;
		double sprt_delta_;
		double sprt_epsilon_;

		// LOSAC parameters
		unsigned int lo_inner_sample_size;
		unsigned int lo_num_inner_ransac_reps_;
		double       lo_threshold_multiplier_;
		unsigned int lo_num_iterative_steps_;

	// ------------------------------------------------------------------------
	// data-specific parameters/temporary storage
	protected:
		// common
		std::vector<unsigned int>     min_sample_;
		std::vector<double>           errs_;
		std::vector<double>::iterator err_ptr_[2];	// err_ptr_[0] points to the temporary error value storage
													// err_ptr_[1] points to the error values for the best solution

		// PROSAC
		unsigned int			  subset_size_prosac_;
		unsigned int			  largest_size_prosac_;
		unsigned int			  stop_len_prosac_;
		std::vector<unsigned int> growth_function_prosac_;
		std::vector<unsigned int> non_random_inliers_prosac_;
		std::vector<unsigned int> maximality_samples_prosac_;

		// SPRT
		double					  decision_threshold_sprt_;
		std::vector<unsigned int> evaluation_pool_;	 // randomized evaluation - holds ordering of points 
		unsigned int			  eval_pool_index_;	 // index to the first point to be verified
		struct TestHistorySPRT
		{
			double epsilon_, delta_, A_;
			unsigned int k_;
			struct TestHistorySPRT *prev_;
		};
		unsigned int last_wald_history_update_;							
		TestHistorySPRT* wald_test_history_;

		// LO
		unsigned int num_prev_best_inliers_lo_;
		std::vector<unsigned int> prev_best_inliers_lo_;

	// ------------------------------------------------------------------------
	// helper functions
	protected:
		// random sampling 
		inline void generateUniformRandomSample(unsigned int dataSize, unsigned int sampleSize, std::vector<unsigned int>* sample);
		inline void generatePROSACMinSample(unsigned int hypCount, std::vector<unsigned int>* sample);

		// PROSAC initialization
		inline void initPROSAC();

		// SPRT updates
		inline void designSPRTTest();
		double computeExpSPRT(double new_epsilon, double epsilon, double delta);
		TestHistorySPRT* addTestHistorySPRT(double epsilon, double delta, unsigned int numHyp, 
											TestHistorySPRT* testHistory, unsigned int* lastUpdate);

		// local optimization
		unsigned int locallyOptimizeSolution(const unsigned int bestInliers);
		unsigned int findInliers(const std::vector<double>::iterator& errPtr, double threshold, 
								 std::vector<unsigned int>* inliers);

		// stopping criterion
		unsigned int updateStandardStopping(unsigned int numInliers, unsigned int totPoints, unsigned int sampleSize);
		unsigned int updatePROSACStopping(unsigned int hypCount);
		unsigned int updateSPRTStopping(unsigned int numInliers, unsigned int totPoints, TestHistorySPRT* testHistory);

		// store USAC result/statistics
		void storeSolution(unsigned int modelIndex, unsigned int numInliers);
};


// ============================================================================================
// initParamsUSAC: initializes USAC 
// this function is called once to initialize the basic USAC parameters 
// all problem specific/data-specific initialization is done using initDataUSAC()
// ============================================================================================
template <class ProblemType>
void USAC<ProblemType>::initParamsUSAC(const ConfigParams& cfg)
{
	// store common parameters
	usac_conf_threshold_        = cfg.common.confThreshold;
	usac_min_sample_size_		= cfg.common.minSampleSize;
	usac_inlier_threshold_		= cfg.common.inlierThreshold;
	usac_max_hypotheses_		= cfg.common.maxHypotheses;
	usac_max_solns_per_sample_	= cfg.common.maxSolutionsPerSample;
	usac_prevalidate_sample_	= cfg.common.prevalidateSample;
	usac_prevalidate_model_		= cfg.common.prevalidateModel;
	usac_test_degeneracy_		= cfg.common.testDegeneracy;
	usac_sampling_method_		= cfg.common.randomSamplingMethod;
	usac_verif_method_			= cfg.common.verifMethod;
	usac_local_opt_method_      = cfg.common.localOptMethod;

	// read in PROSAC parameters if required
	if (usac_sampling_method_ == USACConfig::SAMP_PROSAC)
	{
		prosac_growth_max_samples_  = cfg.prosac.maxSamples;
		prosac_beta_				= cfg.prosac.beta;
		prosac_non_rand_conf_		= cfg.prosac.nonRandConf;
		prosac_min_stop_length_     = cfg.prosac.minStopLen;
	}

	// read in SPRT parameters if required
	if (usac_verif_method_ == USACConfig::VERIF_SPRT)
	{
		sprt_tM_	   = cfg.sprt.tM;
		sprt_mS_       = cfg.sprt.mS;
		sprt_delta_    = cfg.sprt.delta;
		sprt_epsilon_  = cfg.sprt.epsilon;
	}

	// read in LO parameters if required
	if (usac_local_opt_method_ == USACConfig::LO_LOSAC)
	{
		lo_inner_sample_size		   = cfg.losac.innerSampleSize;
		lo_num_inner_ransac_reps_      = cfg.losac.innerRansacRepetitions;
		lo_threshold_multiplier_       = cfg.losac.thresholdMultiplier;
		lo_num_iterative_steps_	       = cfg.losac.numStepsIterative;
	}

	min_sample_.clear(); min_sample_.resize(usac_min_sample_size_);
}

// ============================================================================================
// initDataUSAC: initializes problem specific/data-specific parameters
// this function is called once per USAC run on new data
// ============================================================================================
template <class ProblemType>
void USAC<ProblemType>::initDataUSAC(const ConfigParams& cfg)
{
	// clear output data
	usac_results_.reset();

	// initialize some data specfic stuff
	usac_num_data_points_ = cfg.common.numDataPoints;
	usac_inlier_threshold_ = cfg.common.inlierThreshold*cfg.common.inlierThreshold;

	// set up PROSAC if required
	if (usac_sampling_method_ == USACConfig::SAMP_PROSAC)
	{
		prosac_sorted_point_indices_.assign(cfg.prosac.sortedPointIndices, cfg.prosac.sortedPointIndices+cfg.common.numDataPoints);
		initPROSAC(); 		// init the PROSAC data structures
	}

	// set up SPRT if required
	if (usac_verif_method_ == USACConfig::VERIF_SPRT)
	{
		sprt_tM_	  = cfg.sprt.tM;
		sprt_mS_      = cfg.sprt.mS;
		sprt_delta_   = cfg.sprt.delta;
		sprt_epsilon_ = cfg.sprt.epsilon;

		last_wald_history_update_ = 0;
		wald_test_history_ = NULL;
		designSPRTTest();
	}

	// set up LO if required
	if (usac_local_opt_method_ == USACConfig::LO_LOSAC)
	{
		num_prev_best_inliers_lo_ = 0;
		prev_best_inliers_lo_.clear();
		prev_best_inliers_lo_.resize(usac_num_data_points_, 0);
	}

	// space to store point-to-model errors
	errs_.clear(); errs_.resize(2*usac_num_data_points_, 0);
	err_ptr_[0] = errs_.begin();
	err_ptr_[1] = errs_.begin()+usac_num_data_points_;

	// inititalize evaluation ordering to a random permutation of 0...num_data_points_-1
	eval_pool_index_ = 0;
	evaluation_pool_.clear(); evaluation_pool_.resize(usac_num_data_points_);
	for (unsigned int i = 0; i < usac_num_data_points_; ++i)
	{
		evaluation_pool_[i] = i;
	}
	std::random_shuffle(evaluation_pool_.begin(), evaluation_pool_.end());

	// storage for results
	usac_results_.inlier_flags_.resize(usac_num_data_points_, 0);
	usac_results_.best_sample_.resize(usac_min_sample_size_, 0);
	// if degeneracy testing option is selected
	if (usac_test_degeneracy_)
	{
		usac_results_.degen_inlier_flags_.resize(usac_num_data_points_, 0);
		usac_results_.degen_sample_.resize(usac_num_data_points_, 0);
		usac_results_.degen_sample_.resize(usac_min_sample_size_, 0);
	}
}


// ============================================================================================
// solve: the main USAC function
// computes a solution to the robust estimation problem
// ============================================================================================
template <class ProblemType>
bool USAC<ProblemType>::solve()
{
	unsigned int adaptive_stopping_count = usac_max_hypotheses_;   // initialize with worst case	
	bool update_sprt_stopping = true;	

	// abort if too few data points
	if (usac_num_data_points_ < usac_min_sample_size_ || (usac_sampling_method_ == USACConfig::SAMP_PROSAC && 
		usac_num_data_points_ < prosac_min_stop_length_))
	{
		return false;
	}

	// timing stuff
	double tick = Timer::getTimestampInSeconds();

	// ------------------------------------------------------------------------
	// main USAC loop
	while (usac_results_.hyp_count_ < adaptive_stopping_count && usac_results_.hyp_count_ < usac_max_hypotheses_)
	{
		++usac_results_.hyp_count_;

		// -----------------------------------------
		// 1. generate sample
		switch (usac_sampling_method_)
		{
			case USACConfig::SAMP_UNIFORM:
			{
				generateUniformRandomSample(usac_num_data_points_, usac_min_sample_size_, &min_sample_);
				break;
			}

			case USACConfig::SAMP_UNIFORM_MM:
			{
				generateUniformRandomSample(usac_num_data_points_, usac_min_sample_size_, &min_sample_);
				break;
			}

			case USACConfig::SAMP_PROSAC:
			{
				generatePROSACMinSample(usac_results_.hyp_count_, &min_sample_);
				break;
			} 
		}
		
		// -----------------------------------------
		// 2. validate sample
		if (usac_prevalidate_sample_)
		{
			// pre-validate sample before testing generating model
			bool valid_sample = static_cast<ProblemType *>(this)->validateSample();
			if (!valid_sample)
			{
				++usac_results_.rejected_sample_count_;
				continue;
			}
		}

		// -----------------------------------------
		// 3. generate model(s)
		unsigned int num_solns = static_cast<ProblemType *>
								 (this)->generateMinimalSampleModels();
		usac_results_.model_count_ += num_solns;

		// -----------------------------------------
		// 4-5. validate + evaluate model(s)
		bool update_best = false;
		for (unsigned int i = 0; i < num_solns; ++i)
		{
			if (usac_prevalidate_model_)
			{
				// pre-validate model before testing against data points
				bool valid_model = static_cast<ProblemType *>(this)->validateModel(i);
				if (!valid_model)
				{
					++usac_results_.rejected_model_count_;
					continue;
				}
			}

			// valid model, perform evaluation
			unsigned int inlier_count, num_points_tested;
			bool good = static_cast<ProblemType *>
				        (this)->evaluateModel(i, &inlier_count, &num_points_tested);

			// update based on verification results
			switch (usac_verif_method_)
			{
				case USACConfig::VERIF_STANDARD:
				{
					usac_results_.total_points_verified_ += usac_num_data_points_;
					// check if best so far
					if (inlier_count > usac_results_.best_inlier_count_)
					{
						update_best = true;
						usac_results_.best_inlier_count_ = inlier_count;
						storeSolution(i, usac_results_.best_inlier_count_);    // store result
					}
					break;
				} // end case standard verification
				
				case USACConfig::VERIF_SPRT:
				{
					if (!good)
					{
						// model rejected
						usac_results_.total_points_verified_ += num_points_tested;
						double delta_new = (double)inlier_count/num_points_tested;
						if (delta_new > 0 && abs(sprt_delta_ - delta_new)/sprt_delta_ > 0.1)
						{
							// update parameters
							wald_test_history_ = addTestHistorySPRT(sprt_epsilon_, sprt_delta_, 
								usac_results_.hyp_count_, wald_test_history_, &last_wald_history_update_);
							sprt_delta_ = delta_new;
							designSPRTTest();
						}
					}
					else
					{
						// model accepted
						usac_results_.total_points_verified_ += usac_num_data_points_;
						if (inlier_count > usac_results_.best_inlier_count_)
						{
							// and best so far
							update_best = true;
							usac_results_.best_inlier_count_ = inlier_count;
							wald_test_history_ = addTestHistorySPRT(sprt_epsilon_, sprt_delta_, 
								usac_results_.hyp_count_, wald_test_history_, &last_wald_history_update_);
							sprt_epsilon_ = (double)usac_results_.best_inlier_count_/usac_num_data_points_;
							designSPRTTest();
							update_sprt_stopping = true;
							storeSolution(i, usac_results_.best_inlier_count_);    // store model
						}
					}
				} // end case sprt

			} // end switch verification method
		} // end evaluating all models for one minimal sample

		// -----------------------------------------
		// 6. if new best model, check for degeneracy
		bool degenerate_model = false, upgrade_model = false, upgrade_successful = false;

		if (update_best && usac_test_degeneracy_)
		{
			std::cout << "Testing for degeneracy (" << usac_results_.best_inlier_count_ << ")" << std::endl;
			static_cast<ProblemType *>(this)->testSolutionDegeneracy(&degenerate_model, &upgrade_model);
			if (degenerate_model && upgrade_model)
			{
				// complete model
				unsigned int upgrade_inliers = static_cast<ProblemType *>(this)->upgradeDegenerateModel();
				if (upgrade_inliers > usac_results_.best_inlier_count_)
				{
					usac_results_.best_inlier_count_ = upgrade_inliers;
					upgrade_successful = true;
				}
			}
		}

		// -----------------------------------------
		// 7. perform local optimization if specified
		if (usac_local_opt_method_ == USACConfig::LO_LOSAC && update_best == true)
		{
			std::cout << "(" << usac_results_.hyp_count_ << ") Performing LO. Inlier count before: " << usac_results_.best_inlier_count_;
			unsigned int lo_inlier_count = locallyOptimizeSolution(usac_results_.best_inlier_count_);
			if (lo_inlier_count > usac_results_.best_inlier_count_)
			{
				usac_results_.best_inlier_count_ = lo_inlier_count;
			}
			std::cout << ", inlier count after: " << lo_inlier_count << std::endl;

			if (num_prev_best_inliers_lo_ < usac_results_.best_inlier_count_)
			{
				// backup the old set of inliers for reference in LO
				num_prev_best_inliers_lo_ = usac_results_.best_inlier_count_;
				prev_best_inliers_lo_ = usac_results_.inlier_flags_;
			}
		}

		// -----------------------------------------
		// 8. update the stopping criterion
		if (update_best)
		{
			// update the number of samples required
			if ( usac_sampling_method_ == USACConfig::SAMP_PROSAC && usac_results_.hyp_count_ <= prosac_growth_max_samples_ )
			{
				adaptive_stopping_count = updatePROSACStopping(usac_results_.hyp_count_); 
			}
			else
			{
				adaptive_stopping_count = updateStandardStopping(usac_results_.best_inlier_count_, usac_num_data_points_, usac_min_sample_size_); 
			}
		}
		// update adaptive stopping count to take SPRT test into account
		if (usac_verif_method_ == USACConfig::VERIF_SPRT && usac_sampling_method_ != USACConfig::SAMP_PROSAC)
		{
			if (usac_results_.hyp_count_ >= adaptive_stopping_count && update_sprt_stopping)
			{
				adaptive_stopping_count = updateSPRTStopping(usac_results_.best_inlier_count_, usac_num_data_points_, wald_test_history_);
				update_sprt_stopping = false;
			}
		}

	} // end the main USAC loop

	// ------------------------------------------------------------------------	
	// output statistics
	double tock = Timer::getTimestampInSeconds();
	std::cout << "Number of hypotheses/models: " << usac_results_.hyp_count_ << "/" << usac_results_.model_count_ << std::endl;
	std::cout << "Number of samples rejected by pre-validation: " << usac_results_.rejected_sample_count_ << std::endl;
	std::cout << "Number of models rejected by pre-validation: " << usac_results_.rejected_model_count_ << std::endl;
	std::cout << "Number of verifications per model: " << 
		(double)usac_results_.total_points_verified_/(usac_results_.model_count_-usac_results_.rejected_model_count_) << std::endl;
	std::cout << "Max inliers/total points: " << usac_results_.best_inlier_count_ << "/" << usac_num_data_points_ << std::endl;

	// ------------------------------------------------------------------------
	// timing stuff
	usac_results_.total_runtime_ = tock - tick;
	std::cout << "Time: " << usac_results_.total_runtime_ << std::endl << std::endl;

	// ------------------------------------------------------------------------
	// clean up
	if (usac_verif_method_ == USACConfig::VERIF_SPRT)
	{
		while (wald_test_history_)
		{
			TestHistorySPRT *temp = wald_test_history_->prev_;
			delete wald_test_history_;
			wald_test_history_ = temp;
		}
	}

	return true;
}


// ============================================================================================
// generateUniformRandomSample: generate random sample uniformly distributed between 
// [0...dataSize-1]
// note that the sample vector needs to be properly sized before calling this function
// ============================================================================================
template <class ProblemType>
void USAC<ProblemType>::generateUniformRandomSample(unsigned int dataSize, unsigned int sampleSize, 
													std::vector<unsigned int>* sample)
{
	unsigned int count=0;
	unsigned int index;
	std::vector<unsigned int>::iterator pos;
	pos = sample->begin();
	do {
		index = rand() % dataSize;
		if (find(sample->begin(), pos, index) == pos)
		{
				(*sample)[count] = index;
				++count;
				++pos;
		}
	} while (count < sampleSize);
}


// ============================================================================================
// initPROSAC: initializes PROSAC
// sets up growth function and stopping criterion 
// ============================================================================================
template <class ProblemType> inline
void USAC<ProblemType>::initPROSAC()
{
	// ------------------------------------------------------------------------
	// growth function

	growth_function_prosac_.clear(); growth_function_prosac_.resize(usac_num_data_points_);
	double T_n;
	unsigned int T_n_p = 1; 
	// compute initial value for T_n
	T_n = prosac_growth_max_samples_;
	for (unsigned int i = 0; i < usac_min_sample_size_; ++i)
	{
		T_n *= (double)(usac_min_sample_size_-i)/(usac_num_data_points_-i);
	}
	// compute values using recurrent relation
	for (unsigned int i = 0; i < usac_num_data_points_; ++i)
	{
		if (i+1 <= usac_min_sample_size_)
		{
			growth_function_prosac_[i] = T_n_p;
			continue;
		}
		double temp = (double)(i+1)*T_n/(i+1-usac_min_sample_size_);
		growth_function_prosac_[i] = T_n_p + (unsigned int)ceil(temp - T_n);
		T_n = temp;
		T_n_p = growth_function_prosac_[i];
	}

	// ------------------------------------------------------------------------
	// initialize the data structures that determine stopping

	// non-randomness constraint
	// i-th entry - inlier counts for termination up to i-th point (term length = i+1)
	non_random_inliers_prosac_.clear();
	non_random_inliers_prosac_.resize(usac_num_data_points_, 0);  
	double pn_i = 1.0;    // prob(i inliers) with subset size n
	for (unsigned int n = usac_min_sample_size_+1; n <= usac_num_data_points_; ++n)
	{
		if (n-1 > 1000)
		{
			non_random_inliers_prosac_[n-1] = non_random_inliers_prosac_[n-2];
			continue;
		}

		std::vector<double> pn_i_vec(usac_num_data_points_, 0);
		// initial value for i = m+1 inliers
		pn_i_vec[usac_min_sample_size_] = (prosac_beta_)*std::pow((double)1-prosac_beta_, (double)n-usac_min_sample_size_-1)*(n-usac_min_sample_size_);
		pn_i = pn_i_vec[usac_min_sample_size_];
		for (unsigned int i = usac_min_sample_size_+2; i <= n; ++i)
		{
			// use recurrent relation to fill in remaining values
			if (i == n)
			{
				pn_i_vec[n-1] = std::pow((double)prosac_beta_, (double)n-usac_min_sample_size_);
				break;
			}
			pn_i_vec[i-1] = pn_i * ((prosac_beta_)/(1-prosac_beta_)) * ((double)(n-i)/(i-usac_min_sample_size_+1));
			pn_i = pn_i_vec[i-1];
		}
		// find minimum number of inliers satisfying the non-randomness constraint
		double acc = 0.0;
		unsigned int i_min = 0;
		for (unsigned int i = n; i >= usac_min_sample_size_+1; --i)
		{
			acc += pn_i_vec[i-1];
			if (acc < 1-prosac_non_rand_conf_)
			{
				i_min = i;
			}
			else
			{
				break;
			}
		}
		non_random_inliers_prosac_[n-1] = i_min;
	}

	// maximality constraint
	// i-th entry - number of samples for pool [0...i] (pool length = i+1)
	maximality_samples_prosac_.clear();
	maximality_samples_prosac_.resize(usac_num_data_points_);	
	for (size_t i = 0; i < usac_num_data_points_; ++i)
	{
		maximality_samples_prosac_[i] = usac_max_hypotheses_;
	}

	// other initializations
	largest_size_prosac_ = usac_min_sample_size_;       // largest set sampled in PROSAC
	subset_size_prosac_ = usac_min_sample_size_;		// size of the current sampling pool
	stop_len_prosac_    = usac_num_data_points_;		// current stopping length
}


// ============================================================================================
// generatePROSACMinSample: generate a minimal sample using ordering information
// uses growth function to determine size of current sampling pool
// ============================================================================================
template <class ProblemType>
void USAC<ProblemType>::generatePROSACMinSample(unsigned int hypCount, std::vector<unsigned int>* sample)
{
	// revert to RANSAC-style sampling if maximum number of PROSAC samples have been tested
	if (hypCount > prosac_growth_max_samples_)
	{
		generateUniformRandomSample(usac_num_data_points_, usac_min_sample_size_, sample);	
		return;
	}

	// if current stopping length is less than size of current pool, use only points up to the stopping length
	if (subset_size_prosac_ > stop_len_prosac_)
	{
		generateUniformRandomSample(stop_len_prosac_, usac_min_sample_size_, sample);	
	}

	// increment the size of the sampling pool if required
	if (hypCount > growth_function_prosac_[subset_size_prosac_-1])
	{
		++subset_size_prosac_;
		if (subset_size_prosac_ > usac_num_data_points_)
		{
			subset_size_prosac_ = usac_num_data_points_;
		}
		if (largest_size_prosac_ < subset_size_prosac_)
		{
			largest_size_prosac_ = subset_size_prosac_;
		}
	}

	// generate PROSAC sample
	generateUniformRandomSample(subset_size_prosac_-1, usac_min_sample_size_-1, sample);
	(*sample)[usac_min_sample_size_-1] = subset_size_prosac_-1;
	for (unsigned int i = 0; i < sample->size(); ++i)
	{
		(*sample)[i] = prosac_sorted_point_indices_[(*sample)[i]];
	}
}


// ============================================================================================
// updatePROSACStopping: determines number of samples required for termination
// checks the non-randomness and maximality constraints
// ============================================================================================
template <class ProblemType> inline
unsigned int USAC<ProblemType>::updatePROSACStopping(unsigned int hypCount)
{
	unsigned int max_samples = maximality_samples_prosac_[stop_len_prosac_-1];
	
	// go through sorted points and track inlier counts
	unsigned int inlier_count = 0;

	// just accumulate the count for the first prosac_min_stop_length_ points
	for (unsigned int i = 0; i < prosac_min_stop_length_; ++i)
	{
		inlier_count += usac_results_.inlier_flags_[prosac_sorted_point_indices_[i]];
	}	

	// after this initial subset, try to update the stopping length if possible
	for (unsigned int i = prosac_min_stop_length_; i < usac_num_data_points_; ++i)
	{
		inlier_count += usac_results_.inlier_flags_[prosac_sorted_point_indices_[i]];

		if (non_random_inliers_prosac_[i] < inlier_count)
		{
			non_random_inliers_prosac_[i] = inlier_count;	// update the best inliers for the the subset [0...i]

			// update the number of samples based on this inlier count
			if ( (i == usac_num_data_points_-1) || 
				 (usac_results_.inlier_flags_[prosac_sorted_point_indices_[i]] && !usac_results_.inlier_flags_[prosac_sorted_point_indices_[i+1]]) )
			{
				unsigned int new_samples = updateStandardStopping(inlier_count, i+1, usac_min_sample_size_);
				if (i+1 < largest_size_prosac_)
				{
					// correct for number of samples that have points in [i+1, largest_size_prosac_-1]
					new_samples += hypCount - growth_function_prosac_[i];
				}

				if (new_samples < maximality_samples_prosac_[i])
				{
					// if number of samples can be lowered, store values and update stopping length
					maximality_samples_prosac_[i] = new_samples;
					if ( (new_samples < max_samples) || ( (new_samples == max_samples) && (i+1 >= stop_len_prosac_) ) )
					{
						stop_len_prosac_ = i+1;
						max_samples = new_samples;
					}
				}
			}
		}
	}
	return max_samples;
}


// ============================================================================================
// designSPRTTest: designs a new SPRT test (i.e., updates the SPRT decision threshold) based 
// on current values of delta and epsilon
// ============================================================================================
template <class ProblemType> inline
void USAC<ProblemType>::designSPRTTest()
{
	double An_1, An, C, K;

	C = (1 - sprt_delta_)*log( (1 - sprt_delta_)/(1-sprt_epsilon_) ) 
		+ sprt_delta_*(log( sprt_delta_/sprt_epsilon_ ));
	K = (sprt_tM_*C)/sprt_mS_ + 1;
	An_1 = K;

	// compute A using a recursive relation
	// A* = lim(n->inf)(An), the series typically converges within 4 iterations
	for (unsigned int i = 0; i < 10; ++i)
    {
		An = K + log(An_1);
		if (An - An_1 < 1.5e-8) 
		{
			break;
		}
	    An_1 = An;
    }

	decision_threshold_sprt_ = An;
}

// ============================================================================================
// locallyOptimizeSolution: iteratively refines the current model based on its inlier set
// ============================================================================================
template <class ProblemType> inline
unsigned int USAC<ProblemType>::locallyOptimizeSolution(const unsigned int bestInliers)
{
	// return if insufficient number of points
	if (bestInliers < 2*lo_inner_sample_size) 
	{
		return 0;
	}

	unsigned int lo_sample_size = std::min(lo_inner_sample_size, bestInliers/2);
	std::vector<unsigned int> sample(lo_sample_size);
	std::vector<unsigned int> orig_inliers(usac_num_data_points_);
	std::vector<unsigned int> iter_inliers(usac_num_data_points_);
	unsigned int num_points_tested;

	// find all inliers less than threshold 
	unsigned int lo_inliers = bestInliers;
	unsigned int temp_inliers = 0;
	findInliers(err_ptr_[1], usac_inlier_threshold_, &orig_inliers);	
#if 0
	// check if there is substantial overlap between the current and the best inlier sets
	// if yes, the local refinement is unlikely to help 
	std::vector<unsigned int> ind_best, ind_curr(bestInliers);
	for (size_t i = 0; i < usac_num_data_points_; ++i)
	{
		if (m_prevBestInliers[i]) ind_best.push_back(i);
	}
	for (size_t i = 0; i < bestInliers; ++i)
	{
		ind_curr[i] = orig_inliers[i];
	}
	std::vector<unsigned int> temp_intersection(bestInliers);
	std::vector<unsigned int>::iterator it = std::set_intersection (ind_best.begin(), ind_best.end(), 
																	ind_curr.begin(), ind_curr.end(), 
																	temp_intersection.begin());
	unsigned int num_elements_intersection = it - temp_intersection.begin();
	std::cout << " (" << num_elements_intersection << ") " ;
	if ((double)num_elements_intersection/bestInliers > 0.95)
	{
		std::cout << "****";
		return 0;
	}
#endif
	++usac_results_.num_local_optimizations_;

	double *weights = new double[usac_num_data_points_];	
	double threshold_step_size = (lo_threshold_multiplier_*usac_inlier_threshold_ - usac_inlier_threshold_)
								  /lo_num_iterative_steps_;	
	// perform number of inner RANSAC repetitions
	for (unsigned int i = 0; i < lo_num_inner_ransac_reps_; ++i)
	{
		// generate non-minimal sample model and find inliers 
		generateUniformRandomSample(bestInliers, lo_sample_size, &sample);
		for (unsigned int j = 0; j < lo_sample_size; ++j)
		{
			sample[j] = orig_inliers[sample[j]];    // we want points only from the current inlier set
		}
		if ( !static_cast<ProblemType *>(this)->generateRefinedModel(sample, lo_sample_size) )
		{
			continue;
		}

		if (! static_cast<ProblemType *>(this)->evaluateModel(0, &temp_inliers, &num_points_tested) )
		{
			continue;
		}
		temp_inliers = findInliers(err_ptr_[0], lo_threshold_multiplier_*usac_inlier_threshold_, &iter_inliers);

		// generate least squares model from all inliers
		if (! static_cast<ProblemType *>(this)->generateRefinedModel(iter_inliers, temp_inliers) )
		{
			continue;
		}

		// iterative (reweighted) refinement - reduce threshold in steps, find new inliers and refit fundamental matrix
		// using weighted least-squares
		for (unsigned int j = 0; j < lo_num_iterative_steps_; ++j)
		{
			if (! static_cast<ProblemType *>(this)->evaluateModel(0, &temp_inliers, &num_points_tested) )
			{
				continue;
			}
			findInliers(err_ptr_[0], (lo_threshold_multiplier_*usac_inlier_threshold_) - (j+1)*threshold_step_size, &iter_inliers);		
			static_cast<ProblemType *>(this)->findWeights(0, iter_inliers, temp_inliers, weights);
			if (! static_cast<ProblemType *>(this)->generateRefinedModel(iter_inliers, temp_inliers, true, weights) )
			{
				continue;
			}
		}

		// find final set of inliers for this round
		if (! static_cast<ProblemType *>(this)->evaluateModel(0, &temp_inliers, &num_points_tested) )
		{
			continue;
		}
		//findInliers(err_ptr_[0], iter_inliers, usac_inlier_threshold_);	

		if (temp_inliers > lo_inliers)
		{
			// store model
			lo_inliers = temp_inliers;
			static_cast<ProblemType *>(this)->storeSolution(0, lo_inliers);
		}
	}

	delete[] weights;
	return lo_inliers;
}	


// ============================================================================================
// findInliers: given an error vector and threshold, returns the indices of inliers
// ============================================================================================
template <class ProblemType> inline
unsigned int USAC<ProblemType>::findInliers(const std::vector<double>::iterator& errPtr, double threshold, 
											std::vector<unsigned int>* inliers)
{
	unsigned int inlier_count = 0;
	for (unsigned int i = 0; i < usac_num_data_points_; ++i)
	{
		if (*(errPtr+i) < threshold)
		{
			(*inliers)[inlier_count] = i;
			++inlier_count;
		}
	}
	return inlier_count;
}


// ============================================================================================
// updateStandardStopping: updates the standard RANSAC stopping criterion
// ============================================================================================
template <class ProblemType> inline
unsigned int USAC<ProblemType>::updateStandardStopping(unsigned int numInliers, unsigned int totPoints,
													   unsigned int sampleSize)
{
	double n_inliers = 1.0;
	double n_pts = 1.0;

	for (unsigned int i = 0; i < sampleSize; ++i)
	{
		n_inliers *= numInliers - i;
		n_pts *= totPoints - i;
	}
	double prob_good_model = n_inliers/n_pts;

	if ( prob_good_model < std::numeric_limits<double>::epsilon() )
	{
		return usac_max_hypotheses_;
	}
	else if ( 1 - prob_good_model < std::numeric_limits<double>::epsilon() )
	{
		return 1;
	}
	else 
	{
		double nusample_s = log(1-usac_conf_threshold_)/log(1-prob_good_model);
		return (unsigned int) ceil(nusample_s);
	}
}


// ============================================================================================
// updateSPRTStopping: updates the stopping criterion accounting for erroneous rejections in 
// the SPRT test
// ============================================================================================
template <class ProblemType> inline
unsigned int USAC<ProblemType>::updateSPRTStopping(unsigned int numInliers, unsigned int totPoints, TestHistorySPRT* testHistory)
{
	double n_inliers = 1.0;
	double n_pts = 1.0;
	double h = 0.0, k = 0.0, prob_reject_good_model = 0.0, log_eta = 0.0;
	double new_eps = (double)numInliers/totPoints;
	TestHistorySPRT* current_test = testHistory;

	for (unsigned int i = 0; i < usac_min_sample_size_; ++i)
	{
		n_inliers *= numInliers - i;
		n_pts *= totPoints - i;
	}
	double prob_good_model = n_inliers/n_pts;

	if ( prob_good_model < std::numeric_limits<double>::epsilon() )
	{
		return usac_max_hypotheses_;
	}
	else if ( 1 - prob_good_model < std::numeric_limits<double>::epsilon() )
	{
		return 1;
	}

	while (current_test != NULL)
	{
		k += current_test->k_;
		h = computeExpSPRT(new_eps, current_test->epsilon_, current_test->delta_);
		prob_reject_good_model = 1/(exp( h*log(current_test->A_) ));
		log_eta += (double) current_test->k_ * log( 1 - prob_good_model*(1-prob_reject_good_model) );
		current_test = current_test->prev_;
	}

	double nusample_s = k + ( log(1-usac_conf_threshold_) - log_eta ) / log( 1-prob_good_model * (1-(1/decision_threshold_sprt_)) );
	return (unsigned int) ceil(nusample_s);	
}


// ============================================================================================
// computeExpSPRT: computes the value of the exponent h_i required to determine the number
// of iterations when using SPRT
// ============================================================================================
template <class ProblemType> inline
double USAC<ProblemType>::computeExpSPRT(double newEpsilon, double epsilon, double delta)
{
	double al, be, x0, x1, v0, v1, h;

	al = log(delta/epsilon);
	be = log( (1-delta)/(1-epsilon) );

	x0 = log( 1/(1-newEpsilon) )/be;
	v0 = newEpsilon * exp(x0 *al);
	x1 = log( (1-2*v0) / (1-newEpsilon) )/be;
	v1 = newEpsilon * exp(x1 * al) + (1-newEpsilon) * exp(x1 * be);
	h = x0 - (x0 - x1)/(1+v0 - v1)*v0;
	return h;
}


// ============================================================================================
// addTestHistorySPRT: store statistics for each SPRT test (since each test is adjusted to 
// reflect current estimates of these parameters)
// this is required to compute the number of samples for termination
// ============================================================================================
template <class ProblemType> inline
typename USAC<ProblemType>::TestHistorySPRT* USAC<ProblemType>::addTestHistorySPRT(double epsilon, double delta, unsigned int numHyp, 
																				   TestHistorySPRT* testHistory, 
																				   unsigned int* lastUpdate)
{
	TestHistorySPRT *new_test_history = new TestHistorySPRT;
	new_test_history->epsilon_ = epsilon;
	new_test_history->delta_ = delta;
	new_test_history->A_ = decision_threshold_sprt_;
	new_test_history->k_ = numHyp - *lastUpdate;
	new_test_history->prev_ = testHistory;
	*lastUpdate = numHyp;

	return new_test_history;
}


// ============================================================================================
// storeSolution: stores inlier flag vector, best sample indices and best inlier count
// also swaps error pointers so that err_ptr_[1] points to the error values for the best model
// calls a model-specific function to store model data as appropriate
// ============================================================================================
template <class ProblemType> inline
void USAC<ProblemType>::storeSolution(unsigned int modelIndex, unsigned int numInliers)
{
	// save the current best set of inliers
	usac_results_.best_inlier_count_ = numInliers;
	std::vector<double>::iterator current_err_array = err_ptr_[0];
	for (unsigned int i = 0; i < usac_num_data_points_; ++i)
	{
		if (*(current_err_array+i) < usac_inlier_threshold_)
			usac_results_.inlier_flags_[i] = 1;
		else
			usac_results_.inlier_flags_[i] = 0;
	}

	// save the current best sample indices
	usac_results_.best_sample_ = min_sample_;

	// and switch the error pointers
	std::swap(err_ptr_[0], err_ptr_[1]);

	// store model
	static_cast<ProblemType *>(this)->storeModel(modelIndex, numInliers);
}

#endif   
