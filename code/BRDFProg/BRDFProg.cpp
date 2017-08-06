//************************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: BRDFProg.cpp
//**
//** Thomas A. Germer
//** Optical Technology Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//** Version: 7.00 (January 2015)
//**
//************************************************************************************
#include "brdf.h"

using namespace std;        // Use unqualified names for Standard C++ library
using namespace SCATMECH;   // Use unqualified names for SCATMECH library

int main(int argv,char** argc)
{
    try {

        // Query user for scattering angles and ranges...
        double thetai = AskUser("Incident Angle ",45.)*deg;

        double thetas_1 = AskUser("Initial Scattering Angle ",45.)*deg;
        double thetas_2 = AskUser("Final Scattering Angle ",45.)*deg;
        double thetas_3 = AskUser("Step Scattering Angle ",1.)*deg;

        double phis_1 = AskUser("Initial Azimuthal Angle ",0.)*deg;
        double phis_2 = AskUser("Final Azimuthal Angle ",180.)*deg;
        double phis_3 = AskUser("Step Azimuthal Angle ",2.)*deg;

        // Get an instance of BRDF_Model...
        BRDF_Model_Ptr model = Get_Model_Ptr();

        // Query user for model parameters...
        model->AskUser();

        // Loop through scattering geometries...
        for (double thetas=thetas_1; thetas<=thetas_2; thetas+=thetas_3) {
            for (double phis=phis_1; phis<=phis_2; phis+=phis_3) {

                // Get the Mueller matrix for scattering...
                MuellerMatrix m=model->Mueller(thetai,thetas,phis,0.);

                // Calculate the Stokes vector for p-polarized incident light...
                StokesVector s=m*StokesVectorUnitP();

                // Print out various light scattering parameters...
                cout << thetas/deg << tab   // Scattering angle (theta)
                     << phis/deg << tab     // Scattering angle (phi)
                     << s.eta()/deg << tab  // Principal angle of polarization
                     << s.DOLP() << tab     // Degree of linear polarization
                     << s.DOCP() << tab     // Degree of circular polarization
                     << s.DOP() << tab      // Degree of polarization
                     << s.I() << endl;      // The intensity (BRDF)
            }
        }

    } catch (exception& e) {

        cerr << e.what() << endl;

    }

    return 0;
}
