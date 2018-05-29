//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: sizedistribution.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "sizedistribution.h"

namespace SCATMECH {

	double Log_Normal_Distribution::pdf(double d)
	{
		SETUP();

		return 1. / d / sigma / sqrt(2.*pi)*exp(-sqr(log(d) - log(median)) / 2. / sqr(sigma));
	}

	double Normal_Distribution::pdf(double d)
	{
		SETUP();
		return 1. / sqrt(2.*pi) / sigma*exp(-sqr(d - mean) / 2. / sqr(sigma));
	}

	double Weibull_Distribution::pdf(double d)
	{
		SETUP();
		return exponent / scale*pow(d / scale, exponent - 1.)*exp(-pow(d / scale, exponent));
	}

	double Modified_Gamma_Distribution::pdf(double d)
	{
		SETUP();

		double r = d / 2.;

		double norm = pow(Lambda, -(1 + mu) / gamma)*tgamma((1 + mu) / gamma) / gamma;

		return pow(r, mu)*exp(-Lambda*pow(r, gamma))/norm;
	}

	double Bimodal_Distribution::pdf(double d)
	{
		SETUP();

		return fractionA*distributionA->pdf(d) + (1 - fractionA)*distributionB->pdf(d);
	}

	double CC1246E_SurfaceParticleSizeDistribution::surfacedensity(double d)
	{
		SETUP();

		// This distribution is not defined below 1 micrometer.
		if (d < 1.) return 0.;

		/// cdf = pow(10., 0.926*(sqr(log10(cleanliness)) - sqr(log10(d))) - 11.);
		return 8.68589E-12 * exp(-0.434294*slope*sqr(log(d)) + 0.402157*sqr(log(cleanliness)))*slope*log(d) / d;
	}

	void Register(const Distribution* x)
	{
		static bool regd = false;
		if (!regd) {
			regd = true;

			Register_Model(Distribution);
			Register_Model(Log_Normal_Distribution);
			Register_Model(Normal_Distribution);
			Register_Model(Modified_Gamma_Distribution);
			Register_Model(Weibull_Distribution);
			Register_Model(Bimodal_Distribution);
			Register_Model(VolumeParticleSizeDistribution);
			Register_Model(Regular_VolumeParticleSizeDistribution);
			Register_Model(SurfaceParticleSizeDistribution);
			Register_Model(Regular_SurfaceParticleSizeDistribution);
			Register_Model(CC1246E_SurfaceParticleSizeDistribution);
		}
	}

	DEFINE_VIRTUAL_MODEL(Distribution, Model, "Base class for diameter distributions");

	DEFINE_VIRTUAL_MODEL(VolumeParticleSizeDistribution, Model, "Volume particle size and number distribution");

	DEFINE_MODEL(Regular_VolumeParticleSizeDistribution, VolumeParticleSizeDistribution, "Volume particle size and number distribution");
	DEFINE_PTRPARAMETER(Regular_VolumeParticleSizeDistribution, Distribution_Ptr, distribution, "Diameter distribution", "Log_Normal_Distribution", 0xFF);
	DEFINE_PARAMETER(Regular_VolumeParticleSizeDistribution, double, numberdensity, "Volume number density [1/um^3]", "0.00001", 0xFF);

	DEFINE_VIRTUAL_MODEL(SurfaceParticleSizeDistribution, Model, "Surface particle size and number distribution");

	DEFINE_MODEL(Regular_SurfaceParticleSizeDistribution, SurfaceParticleSizeDistribution, "Surface particle size and number distribution");
	DEFINE_PTRPARAMETER(Regular_SurfaceParticleSizeDistribution, Distribution_Ptr, distribution, "Diameter distribution", "Log_Normal_Distribution", 0xFF);
	DEFINE_PARAMETER(Regular_SurfaceParticleSizeDistribution, double, numberdensity, "Surface number density [1/um^2]", "0.001", 0xFF);

	DEFINE_MODEL(Log_Normal_Distribution, Distribution, "Log-Normal distribution");
	DEFINE_PARAMETER(Log_Normal_Distribution, double, sigma, "Standard deviation of log(diameter)", "1", 0xFF);
	DEFINE_PARAMETER(Log_Normal_Distribution, double, median, "Median diameter [um]", "10", 0xFF);

	DEFINE_MODEL(Modified_Gamma_Distribution, Distribution, "Modified gamma distribution of Deirmendjian");
	DEFINE_PARAMETER(Modified_Gamma_Distribution, double, mu, "mu", "8", 0xFF);
	DEFINE_PARAMETER(Modified_Gamma_Distribution, double, Lambda, "Lambda", "0.0417", 0xFF);
	DEFINE_PARAMETER(Modified_Gamma_Distribution, double, gamma, "Gamma Parameter", "3", 0xFF);

	DEFINE_MODEL(Weibull_Distribution, Distribution, "Weibull distribution");
	DEFINE_PARAMETER(Weibull_Distribution, double, scale, "Scale (lambda) [um]", "10", 0xFF);
	DEFINE_PARAMETER(Weibull_Distribution, double, exponent, "Exponent (k)", "4", 0xFF);
	
	DEFINE_MODEL(Normal_Distribution, Distribution, "Normal distribution");
	DEFINE_PARAMETER(Normal_Distribution, double, mean, "Mean diameter [um]", "10", 0xFF);
	DEFINE_PARAMETER(Normal_Distribution, double, sigma, "Standard deviation [um]", "5", 0xFF);

	DEFINE_MODEL(Bimodal_Distribution, Distribution, "Sum of two distributions");
	DEFINE_PTRPARAMETER(Bimodal_Distribution, Distribution_Ptr, distributionA, "First distribution", "Log_Normal_Distribution", 0xFF);
	DEFINE_PTRPARAMETER(Bimodal_Distribution, Distribution_Ptr, distributionB, "Second distribution", "Log_Normal_Distribution", 0xFF);
	DEFINE_PARAMETER(Bimodal_Distribution, double, fractionA, "Fraction of first distribution", "0.5", 0xFF);

	DEFINE_MODEL(CC1246E_SurfaceParticleSizeDistribution, SurfaceParticleSizeDistribution,"Size distribution defined by cleanliness standard IEST-STD-CC1246E");
	DEFINE_PARAMETER(CC1246E_SurfaceParticleSizeDistribution, double, cleanliness, "Cleanliness level", "1000", 0xFF);
	DEFINE_PARAMETER(CC1246E_SurfaceParticleSizeDistribution, double, slope, "Slope", "0.926", 0xFF);

}
