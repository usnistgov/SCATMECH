//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: torrspar.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//** Version: 7.00 (January 2015)
//**
//******************************************************************************

#include "scatmech.h"
#include "torrspar.h"
#include "vector3d.h"

using namespace std;


namespace SCATMECH {


    //
    // Scattering matrix for Shadowed_Facet_BRDF_Model...
    //
    JonesMatrix
    Shadowed_Facet_BRDF_Model::
    jones()
    {
        // This program implements the model outlined in
        // K.E.Torrance and E.M.Sparrow, "Theory of off-specular reflection
        // from roughened surfaces,"
        // J. Opt. Soc. Am. vol. 57, no. 9, pp. 1105-1114 (1967).

        // All BRDF_Model::jones(,,) should call setup if recalc is set...
        SETUP();

        double s = shadow->f(thetai,thetas,phis);

        // The scattering matrix is the Facet_BRDF_Model scattering
        // matrix multiplied by the T&S shadow function...
        return Facet_BRDF_Model::jones()*
               (COMPLEX)(sqrt(s));
    }

    static
    double min3(double x,double y, double z)
    {
        if (x<y) {
            if (x<z) return x;
            else return z;
        } else {
            if (y<z) return y;
            else return z;
        }
    }

    double
    Torrance_Sparrow_Shadow_Function::
    f(double thetai,double thetas,double phis)
    {
        // The article by Torrance and Sparrow gives a more complicated
        // version of the shadowing function.  It is significantly
        // simplified in R.L. Cook and K.E. Torrance, "A Reflectance
        // Model for Computer Graphics,"
        // Computer Graphics vol. 15, no. 3, pp. 307-316 (1981)
        //
        // In SCATMECH Version 3.00, we use the simplified expressions...
        //
        double cos_thetas = cos(thetas);
        double cos_thetai = cos(thetai);
        double sin_thetas = sin(thetas);
        double sin_thetai = sin(thetai);
        double cos_phis = cos(phis);
        double sin_phis = sin(phis);

        Vector in(-sin_thetai,0.,cos_thetai);
        Vector out(sin_thetas*cos_phis,sin_thetas*sin_phis,cos_thetas);
        Vector normal = unit(in+out);
        double cos_delta = normal.z;
        double cos_beta = normal*in;
        double temp=2.*cos_delta/cos_beta;
        return min3(1.,temp*cos_thetai,temp*cos_thetas);
    }

    double
    Table_Shadow_Function::
    f(double thetai,double thetas,double phis)
    {
        return T.value(thetai)*T.value(thetas);
    }

    static double erfc(double x);

    double
    Smith_Shadow_Function::
    f(double thetai,double thetas,double phis)
    {
        SETUP();

        return S(thetai)*S(thetas);
    }

    double
    Smith_Shadow_Function::
    S(double theta)
    {
        if (theta!=0.) {
            const double sqrt2 = 1.4142135623730951;
            const double sqrttowoverpi = 0.7978845608028654;
            double tant = tan(fabs(theta));
            double mu = 1./tant;
            double Lambda = 0.5*(sqrttowoverpi*w/mu*exp(-sqr(mu/w)/2)-erfc(mu/w/sqrt2));
            return (1.-0.5*erfc(mu/w/sqrt2))/(Lambda+1.);
        } else {
            return 1.;
        }
    }


    double
    Maxwell_Beard_Shadow_Function::
    f(double thetai,double thetas,double phis)
    {
        SETUP();

        double temp = sqr(sin(thetai))-2.*sin(thetai)*sin(thetas)*cos(phis)+
                      sqr(sin(thetas));

        double localslope = (temp>0.) ? sqrt(temp)/(cos(thetai)+cos(thetas)) : 0.;

        double thetan = atan(localslope);

        double cos_angle = sqrt((1.-sin(thetai)*sin(thetas)*cos(phis)+
                                 cos(thetai)*cos(thetas))/2.);
        double beta = (cos_angle<1.) ? acos(cos_angle) : 0.;

        return (1.+(thetan/Omega)*exp(-2*beta/tau))/(1.+(thetan/Omega));
    }


    static double erf(double xx)
    {
        double x = fabs(xx);
        double y;
        double x2 = x*x, x3 = x2*x, x4 = x3*x;

        if (x<0.8) {
            if (x<0.4) {
                if (x<0.2) {
                    if (x<0.1) {
                        // 0.0 < x < 0.1
                        y = 3.9991158574617776e-9 + 1.1283778767060675*x + 0.00009222917635105451*x2 -
                            0.3786078939817188*x3 + 0.0280596900835568*x4;
                    } else {
                        // 0.1 < x < 0.2
                        y = 6.9123302247420595e-6 + 1.1281348191596938*x + 0.0034155717156469922*x2 -
                            0.3997714699494086*x3 + 0.08141147871092347*x4;
                    }
                } else {
                    if (x<0.3) {
                        // 0.2 < x < 0.3
                        y = 0.00008289807803120094 + 1.1266611631809869*x + 0.014269800057487458*x2 -
                            0.43578435264256754*x3 + 0.12684013745036538*x4;
                    } else {
                        // 0.3 < x < 0.4
                        y = 0.0003515940764872161 + 1.1231042881517435*x + 0.03203302755250448*x2 -
                            0.4754525263585343*x3 + 0.16026254442919985*x4;
                    }
                }
            } else {
                if (x<0.6) {
                    if (x<0.5) {
                        // 0.4 < x < 0.5
                        y = 0.0008142096889967032 + 1.1184485911289581*x + 0.04966280414399904*x2 -
                            0.5052223160451739*x3 + 0.17917536429668068*x4;
                    } else {
                        // 0.5 < x < 0.6
                        y = 0.0009894040356612788 + 1.1169378347505283*x + 0.05453332148339385*x2 -
                            0.5121827196508875*x3 + 0.18289705227883957*x4;
                    }
                } else {
                    if (x<0.7) {
                        // 0.6 < x < 0.7
                        y = -0.0004590679411458076 + 1.126401065645961*x + 0.03132313121263053*x2 -
                            0.48685288113118474*x3 + 0.17251862138838447*x4;
                    } else {
                        // 0.7 < x < 0.8
                        y = -0.005868435225691471 + 1.1570536164370822*x - 0.033880559439282365*x2 -
                            0.4251427453415104*x3 + 0.150593633724025*x4;
                    }
                }
            }
        } else {
            if (x<1.5) {
                if (x<1) {
                    if (x<0.9) {
                        // 0.8 < x < 0.9
                        y = -0.018299612016534184 + 1.2189259886210309*x - 0.1494584793296525*x2 -
                            0.3291063473599394*x3 + 0.12064370142014269*x4;
                    } else {
                        // 0.9 < x < 1.0
                        y = -0.040789152550122054 + 1.3186213975822056*x - 0.3152994045856019*x2 -
                            0.20641389445896152*x3 + 0.08658184946095915*x4;
                    }
                } else {
                    if (x<1.25) {
                        // 1.0 < x < 1.25
                        y = 0.2939123365159174 + 1.2433693310809057*x - 0.28005534093283835*x2 +
                            0.12416185515298261*x3 - 0.02138609756883758*x4;
                    } else {
                        // 1.25 < x < 1.5
                        y = 0.32161653107810695 + 1.1539262516648403*x - 0.1715260962994689*x2 +
                            0.06550539608740369*x3 - 0.009472397344441497*x4;
                    }
                }
            } else {
                if (x<3) {
                    if (x<2) {
                        // 1.5 < x < 2.0
                        y = 0.35500102962660396 + 1.0669130364930557*x - 0.08621980792835265*x2 +
                            0.02822007954874106*x3 - 0.003342151491744008*x4;
                    } else {
                        // 2.0 < x < 3.0
                        y = 0.4017842025285597 + 0.9752238499006763*x - 0.01838486624730159*x2 +
                            0.005767242280957822*x3 - 0.0005371836666146841*x4;
                    }
                } else {
                    if (x<6) {
                        // 3.0 < x < 0.6
                        y = 0.08701217120068688 + 1.2660535332655487*x - 0.11460752190343908*x2 +
                            0.018909030317161776*x3 - 0.0011198336995586978*x4;
                    } else {
                        // 6 < x
                        if (xx<0) return -1.;
                        else return 1.;
                    }
                }
            }
        }
        if (x>1) y = 1-exp(-y*y);
        if (xx<0) return -y;
        else return y;
    }

    static double erfc(double x) {
        return 1.-erf(x);
    }


    void Register(const Shadow_Function* x)
    {
        static bool Models_Registered = false;
        if (!Models_Registered) {
            Models_Registered=true;

            Register_Model(Shadow_Function);
            Register_Model(Unit_Shadow_Function);
            Register_Model(Torrance_Sparrow_Shadow_Function);
            Register_Model(Smith_Shadow_Function);
            Register_Model(Maxwell_Beard_Shadow_Function);
            Register_Model(Table_Shadow_Function);
        }
    }

    DEFINE_VIRTUAL_MODEL(Shadow_Function,Model,
                         "Shadow Function");

    DEFINE_MODEL(Unit_Shadow_Function,Shadow_Function,
                 "Unit Shadow Function");

    DEFINE_MODEL(Torrance_Sparrow_Shadow_Function,Shadow_Function,
                 "Shadow function appropriate for V-grooves");

    DEFINE_MODEL(Shadowed_Facet_BRDF_Model,Facet_BRDF_Model,
                 "The facet scattering model with a shadow function.");

    DEFINE_MODEL(Smith_Shadow_Function,Shadow_Function,
                 "Shadow function developed by Smith appropriate for gaussian rough surfaces");
    DEFINE_PARAMETER(Smith_Shadow_Function,double,w,"RMS surface slope","0.35",0xFF);

    DEFINE_MODEL(Maxwell_Beard_Shadow_Function,Shadow_Function,
                 "Empirical shadow function developed by Maxwell-Beard");
    DEFINE_PARAMETER(Maxwell_Beard_Shadow_Function,double,Omega,"Omega [rad]","0.7",0xFF);
    DEFINE_PARAMETER(Maxwell_Beard_Shadow_Function,double,tau,"tau [rad]","0.25",0xFF);

    DEFINE_MODEL(Table_Shadow_Function,Shadow_Function,
                 "Shadow function from a table.");

    DEFINE_PTRPARAMETER(Shadowed_Facet_BRDF_Model,Shadow_Function_Ptr,shadow,"Shadow function","Torrance_Sparrow_Shadow_Function",0xFF);

    DEFINE_PARAMETER(Table_Shadow_Function,Table,T,"Tabulated shadow function","1",0xFF);


} // namespace SCATMECH

