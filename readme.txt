USAC Version 1.0
=================

Overview:
--------

USAC, or Universal RANSAC, is a framework for robust estimation. USAC extends the simple hypothesize-and-verify structure of standard RANSAC to incorporate a number of important practical and computational considerations. This implementation thus addresses many of the limitations of standard RANSAC (efficiency, robustness, degeneracy, accuracy) within a single unified package.

A more in-depth discussion may be found in the following paper:

R. Raguram, O. Chum, M. Pollefeys, J. Matas, J-M. Frahm. "USAC: A Universal Framework for Random Sample Consensus". PAMI (under review)

------------------
Using the library
------------------

This initial release contains a version of USAC that has been tested on Windows, with Visual Studio 2008/2010. To run the code, build the solution file located under msvc/USAC. The 'RunSingleTest' application provides a simple introduction to using the library. 

In brief, calling the estimator consists of these basic steps (shown below for homography estimation):
		
		HomogEstimator* homog = new HomogEstimator;

		// initialize the USAC parameters, either from a config file, or from your application
		homog->initParamsUSAC(cfg);   

		// get input data points/prosac ordering data (if required)
		// set up point_data, cfg.common.numDataPoints, cfg.prosac.sortedPointIndices
		
		// set up the estimation problem
		homog->initDataUSAC(cfg);
		homog->initProblem(cfg, &point_data);
		
		// solve
		if (!homog->solve())
		{
			return(EXIT_FAILURE);
		}

		// do stuff
		
		// cleanup
		homog->cleanupProblem();
		delete homog;

Configuring USAC:
-----------------
The data/ folder contains example configuration files that may be used to control the operation of USAC. Individual modules (e.g., non-uniform sampling with PROSAC) may be switched on or off in this file. Turning on all optimizations provides the best performance, corresponding to USAC-1.0. The library may also be configured to operate as a specific RANSAC variant (e.g., LO-RANSAC), by turning on only the appropriate option in the config file. 

Dependencies (ext/):
--------------------
USAC currently requires the following external libraries to be provided:
LAPACK: For linear algebra routines (used in the fundamental matrix/homography estimators)
LIBCONFIG: Used to configure the operation of USAC (used to demonstrate the usage of parameters)

Solving a specific estimation problem with USAC:
-------------------------------------------------
The src/estimators folder contains the main USAC estimator, along with sample implementations of fundamental matrix and homography estimators. Thus, the main robust estimator is decoupled from the specific estimation problem. 

To write your own estimator to solve a different problem (e.g., fitting planes to 3D points), you need to derive a class from the base USAC class, and implement the data-specific functionality (e.g., computing model paramters from minimal samples, etc.) within it. In particular, there should be an implementation for each of the pure virtual functions in USAC.h (however, many of these can be empty). 

To help with this process, a template (TemplateEstimator.h) is provided in the src/estimators folder, which can be used as a basis for new estimation problems. In most cases, you should not have to change anything in USAC.h, unless you want to alter the behaviour of the robust estimation algorithm itself. 

Misc:
-----
The matlab/ folder has some simple tools that can be used to visualize the output of the estimation process:
    * show_inliers: reads in an image pair, and draws inlier correspondences between the images
	* show_h_transformed: for homography estimation; warps image into alignment based on the estimated homography
	* show_epipolar_lines: for fundamental matrix estimation; steps through epipolar lines corresponding to the computed epipolar geometry 

Bug reports:
------------
Please contact rraguram@cs.unc.edu with bug reports, or questions about this library. 