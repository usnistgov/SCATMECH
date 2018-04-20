//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: firstdiffuse.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "firstdiffuse.h"

#include "askuser.h"
#include "matrix3d.h"

using namespace std;


namespace SCATMECH {


    //
    // Constructor...
    //
    First_Diffuse_BRDF_Model::
    First_Diffuse_BRDF_Model()
    {
        l_scat = 1.;
        phase_function = new Isotropic_Phase_Function;
    }

    void
    First_Diffuse_BRDF_Model::
    setup()
    {
        BRDF_Model::setup();
    }

    JonesMatrix
    First_Diffuse_BRDF_Model::
    jones()
    {
        SETUP();

        throw_backward();
        throw_transmission();

        if (is_reflection()) {
            double n = substrate.n(lambda);

            // Use the absorption coefficient (k) to determine the absorption rate...
            double k_absorb = 4*pi/lambda*substrate.k(lambda);
            double atten = (k_absorb)/(1./l_scat+k_absorb);

            // Angles of refraction internal to the material...
            double sin_thetai_internal = sin(thetai)/n;
            double cos_thetai_internal = sqrt(1.-sqr(sin_thetai_internal));
            double sin_thetas_internal = sin(thetas)/n;
            double cos_thetas_internal = sqrt(1.-sqr(sin_thetas_internal));

            // The incident and scattering khat vectors outside the material...
            Vector inKo(sin(thetai),0.,-cos(thetai));
            Vector outKo(sin(thetas)*cos(phis),sin(thetas)*sin(phis),cos(thetas));

            // The incident and scattering khat vectors inside the material...
            // (one of my compilers not have a COMPLEX sin(COMPLEX)!)
            Vector inKi(sin_thetai_internal,0.,-cos_thetai_internal);
            Vector outKi(sin_thetas_internal*cos(phis),
                         sin_thetas_internal*sin(phis),
                         cos_thetas_internal);

            // The surface normal...
            Vector normal(0.,0.,1.);

            // The following definitions make {s,p,k} right handed...
            Vector cnormal = normal;
            Vector inSi =  perpto(inKi,cnormal);
            Vector inPi =  perpto(inKi,inSi);
            Vector outSi = perpto(outKi,cnormal);
            Vector outPi = perpto(outKi,outSi);
            Vector inSo =  perpto(inKo,cnormal);
            Vector inPo =  perpto(inKo,inSo);
            Vector outSo = perpto(outKo,cnormal);
            Vector outPo = perpto(outKo,outSo);

            // The transmittance for incident light...
            double field_power_i = sqrt(n*cos_thetai_internal/cos(thetai));
            COMPLEX tsi = stack->ts12(thetai,lambda,vacuum,substrate)*field_power_i;
            COMPLEX tpi = stack->tp12(thetai,lambda,vacuum,substrate)*field_power_i;

            // The transmittance for scattered light...
            double field_power_s = sqrt(n*cos_thetas_internal/cos(thetas));
            COMPLEX tss = stack->ts12(thetas,lambda,vacuum,substrate)*field_power_s;
            COMPLEX tps = stack->tp12(thetas,lambda,vacuum,substrate)*field_power_s;

            // Transmission transfer matrices...
            CMatrix ti = (tsi*(CMatrix)outer(inSi,inSo))+(tpi*(CMatrix)outer(inPi,inPo));
            CMatrix ts = (tss*(CMatrix)outer(outSo,outSi))+(tps*(CMatrix)outer(outPo,outPi));

            // The Jacobian converting internal solid angle to external solid angle...
            double jacobian = cos(thetas)/cos_thetas_internal/sqr(n);

            // The cosine of the internal, local scattering angle...
            double cos_scatterangle = inKi*outKi;

            // The probability of scattering that angle per solid angle, if the scattering were scalar...
            // For g=0, the scatter function will be Rayleigh-like...
            double scatter = phase_function->f(acos(cos_scatterangle));

            // Normalize that probability, because Rayleigh scattering is not isotropic...
            //scatter *= 1.5;
            scatter /= 0.5*(1.+sqr(cos_scatterangle));

            // The Rayleigh scattering function...
            CMatrix scatterlocal = (outer(outPi,outPi)+outer(outSi,outSi));

            // The probability that a ray is scattered, and that that scattered ray reaches
            // the surface, is given by...
            double first = sqr(1.-atten)*cos_thetas_internal/(cos_thetai_internal+cos_thetas_internal);

            // The scatter matrix globally includes transmission to/from defect...
            CMatrix scatterglobal = ts*(scatterlocal*ti);

            COMPLEX pp = (CVector)outPo*(scatterglobal*(CVector)inPo);
            COMPLEX sp = (CVector)outPo*(scatterglobal*(CVector)inSo);
            COMPLEX ps = (CVector)outSo*(scatterglobal*(CVector)inPo);
            COMPLEX ss = (CVector)outSo*(scatterglobal*(CVector)inSo);

            // Collect common factors, including the cos(thetas) required for BRDF...
            COMPLEX common=scatter*first*jacobian/cos(thetas);

            // Return the whole value with factors common to all elements...
            return JonesMatrix(pp,ss,ps,sp)*sqrt(common);

        }
        return JonesZero();
    }

    double
    Henyey_Greenstein_Phase_Function::
    f(double theta)
    {
        SETUP();
        return (1.-sqr(g))/pow(1.+sqr(g)-2.*g*cos(theta),1.5)/4./pi;
    }

    double
    Double_Henyey_Greenstein_Phase_Function::
    f(double theta)
    {
        SETUP();
        double q1 = (1.-sqr(g1))/pow(1.+sqr(g1)-2.*g1*cos(theta),1.5)/4./pi;
        double q2 = (1.-sqr(g2))/pow(1.+sqr(g2)-2.*g2*cos(theta),1.5)/4./pi;
        return fabs((1+c)/2.*q1+(1-c)/2.*q2);
    }

    double
    Table_Phase_Function::
    f(double theta)
    {
        SETUP();
        return table.value(theta/deg);
    }

    void
    Reynolds_McCormick_Phase_Function::
    setup()
    {
        Phase_Function::setup();
        double temp1 = pow(1-sqr(g),2*alpha);
        double temp2 = pow(1+sqr(g),2*alpha);
        K = alpha*g/pi*temp1/(temp2-temp1);
    }

    double
    Reynolds_McCormick_Phase_Function::
    f(double theta)
    {
        SETUP();
        return K*pow(1+sqr(g)-2*g*cos(theta),-(alpha+1));
    }

    void
    Double_Reynolds_McCormick_Phase_Function::
    setup()
    {
        Phase_Function::setup();
        double temp1 = pow(1-sqr(g),2*alpha);
        double temp2 = pow(1+sqr(g),2*alpha);
        K = alpha*g/pi*temp1/(temp2-temp1);
    }

    double
    Double_Reynolds_McCormick_Phase_Function::
    f(double theta)
    {
        SETUP();
        double q1 = K*pow(1+sqr(g)-2*g*cos(theta),-(alpha+1));
        double q2 = K*pow(1+sqr(g)+2*g*cos(theta),-(alpha+1));
        return fabs((1+c)/2.*q1+(1-c)/2.*q2);
    }

	static double LegendreP(int i,double x) {
		if (i<0) return 0;
		if (i==0) return 1;
		if (i==1) return x;
		double P1=x;
		double P2=1;
		for (int j=2;j<=i;++j) {
			double P = ((2.*j-1.)*x*P1-(j-1.)*P2)/j;
			P2=P1;
			P1=P;
		}
		return P1;
	}

	double 
	Legendre_Phase_Function::
	f(double theta)
	{
		SETUP();
		double result=0;
		double costheta=cos(theta);
		if (c0!=0) result += c0*LegendreP(0,costheta);
		if (c1!=0) result += 3*c1*LegendreP(1,costheta);
		if (c2!=0) result += 5*c2*LegendreP(2,costheta);
		if (c3!=0) result += 7*c3*LegendreP(3,costheta);
		if (c4!=0) result += 9*c4*LegendreP(4,costheta);
		if (c5!=0) result += 11*c5*LegendreP(5,costheta);
		return result;
	}

	double Gaussian_Phase_Function::f(double theta)
	{
		SETUP();
		double wrad = width*deg;
		double result = (1 - background)*exp(-sqr(theta / wrad) / 2) + background;
		// The following is an approximation good for small widths ...
		double wrad2 = sqr(wrad);
		double wrad4 = sqr(wrad2);
		double wrad6 = wrad4*wrad2;
		double norm = 2 * background + (1 - background) * wrad2 + (background - 1) / 3 * wrad4 + (1 - background) / 15 * wrad6;
		return result / norm / 2 / pi;
	}


    void Register(const First_Diffuse_BRDF_Model* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(First_Diffuse_BRDF_Model);
            Register((Phase_Function*)0);
        }
    }

    void Register(const Phase_Function* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(Phase_Function);
            Register_Model(Henyey_Greenstein_Phase_Function);
            Register_Model(Double_Henyey_Greenstein_Phase_Function);
            Register_Model(Reynolds_McCormick_Phase_Function);
            Register_Model(Double_Reynolds_McCormick_Phase_Function);
            Register_Model(Isotropic_Phase_Function);
            Register_Model(Legendre_Phase_Function);
            Register_Model(Rayleigh_Phase_Function);
			Register_Model(Gaussian_Phase_Function);
			Register_Model(Table_Phase_Function);
        }
    }

    DEFINE_MODEL(First_Diffuse_BRDF_Model,BRDF_Model,
                 "Model for single scattering in a diffuse material under a smooth interface.");

    DEFINE_VIRTUAL_MODEL(Phase_Function,Model,
                         "General phase function for volume scattering");

    DEFINE_MODEL(Henyey_Greenstein_Phase_Function,Phase_Function,
                 "Henyey-Greenstein phase function");

    DEFINE_MODEL(Double_Henyey_Greenstein_Phase_Function,Phase_Function,
                 "Double Henyey-Greenstein phase function with retroreflection");

    DEFINE_MODEL(Isotropic_Phase_Function,Phase_Function,
                 "Scattering which is isotropic.");

    DEFINE_MODEL(Rayleigh_Phase_Function,Phase_Function,
                 "Phase function appropriate for Rayleigh scattering");

    DEFINE_MODEL(Table_Phase_Function,Phase_Function,
                 "Phase function provided by a table of values versus angle");

    DEFINE_MODEL(Reynolds_McCormick_Phase_Function,Phase_Function,
                 "Phase function described by Reynolds and McCormick, J. Opt. Soc. Am. vol. 70, pp 1206-1212 (1980)");

    DEFINE_MODEL(Double_Reynolds_McCormick_Phase_Function,Phase_Function,
                 "A Reynolds-McCormick phase function with a opposition effect.");

    DEFINE_MODEL(Legendre_Phase_Function,Phase_Function,
                 "Phase function defined by a series of up to six Legendre polynomials (l=0 to l=5)");

	DEFINE_MODEL(Gaussian_Phase_Function, Phase_Function,
		"Phase function defined by a Gaussian");

    DEFINE_PARAMETER(First_Diffuse_BRDF_Model,double,l_scat,"Scatter mean free path [um]","1",0xFF);

    DEFINE_PTRPARAMETER(First_Diffuse_BRDF_Model,Phase_Function_Ptr,phase_function,"Phase Function","Henyey_Greenstein_Phase_Function",0xFF);

    DEFINE_PTRPARAMETER(First_Diffuse_BRDF_Model,StackModel_Ptr,stack,"Film stack on substrate","No_StackModel",0xFF);

    DEFINE_PARAMETER(Henyey_Greenstein_Phase_Function,double,g,"Asymmetry parameter (g)","0.01",0xFF);

    DEFINE_PARAMETER(Double_Henyey_Greenstein_Phase_Function,double,g1,"Asymmetry parameter (g1)","0.01",0xFF);
	DEFINE_PARAMETER(Double_Henyey_Greenstein_Phase_Function,double,g2, "Asymmetry parameter (g2)", "0.01", 0xFF);
	DEFINE_PARAMETER(Double_Henyey_Greenstein_Phase_Function,double,c,"Parameter c (-1<c<1)","0.1",0xFF);

    DEFINE_PARAMETER(Table_Phase_Function,Table,table,"Tabulated phase function","0.07957747",0xFF);

    DEFINE_PARAMETER(Reynolds_McCormick_Phase_Function,double,g,"Asymmetry parameter (g)","0.01",0xFF);
    DEFINE_PARAMETER(Reynolds_McCormick_Phase_Function,double,alpha,"Peaking factor (alpha)","0.5",0xFF);

    DEFINE_PARAMETER(Double_Reynolds_McCormick_Phase_Function,double,g,"Asymmetry parameter (g)","0.01",0xFF);
    DEFINE_PARAMETER(Double_Reynolds_McCormick_Phase_Function,double,alpha,"Peaking factor (alpha)","0.5",0xFF);
    DEFINE_PARAMETER(Double_Reynolds_McCormick_Phase_Function,double,c,"Mixing Factor (c)","0.1",0xFF);

    DEFINE_PARAMETER(Legendre_Phase_Function,double,c0,"Legendre coefficient (l=0)","1",0xFF);
    DEFINE_PARAMETER(Legendre_Phase_Function,double,c1,"Legendre coefficient (l=1)","0.65",0xFF);
    DEFINE_PARAMETER(Legendre_Phase_Function,double,c2,"Legendre coefficient (l=2)","0.42",0xFF);
    DEFINE_PARAMETER(Legendre_Phase_Function,double,c3,"Legendre coefficient (l=3)","0",0xFF);
    DEFINE_PARAMETER(Legendre_Phase_Function,double,c4,"Legendre coefficient (l=4)","0",0xFF);
    DEFINE_PARAMETER(Legendre_Phase_Function,double,c5,"Legendre coefficient (l=5)","0",0xFF);

	DEFINE_PARAMETER(Gaussian_Phase_Function, double, width, "Width [deg] parameter", "20", 0xFF);
	DEFINE_PARAMETER(Gaussian_Phase_Function, double, background, "Background parameter", "0.1", 0xFF);

} // namespace SCATMECH
