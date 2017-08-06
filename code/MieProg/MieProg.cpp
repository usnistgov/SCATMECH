//************************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: MieProg.cpp
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
#include "miescat.h"

using namespace std;        // Use unqualified names for Standard C++ library
using namespace SCATMECH;   // Use unqualified names for the SCATMECH library

int main(int argv,char** argc)
{
    try {

        // Query user for angle and wavelength...
        double dtheta = AskUser("Step angle",1.)*deg;

        // Create a Mie scatterer...
        MieScatterer mie;

        // Query user for parameters...
        mie.AskUser();

        // Output scattering efficiencies...
        cout << "Qsca = " << mie.Qsca()
             << "  Qext = " << mie.Qext()
             << "  Qback = " << mie.Qback() << endl;

        // Output column header...
        cout << endl
             << "Angle" << tab << "S11" << tab << "Pol" << tab << "S33" << tab << "S34" << endl;

        // Get normalization factor at zero angle...
        double m00= MuellerMatrix(mie.s(0.0))[0][0];

        // Loop through angles...
        for (double theta=0.; theta<180.*deg; theta+=dtheta) {

            // Evaluate scattering matrix...
            MuellerMatrix m = mie.s(theta);

            // Output unique elements of normalized Mueller matrix...
            cout << theta/deg << tab
                 << m[0][0]/m00 << tab
                 << m[0][1]/m[0][0] << tab
                 << m[2][2]/m[0][0] << tab
                 << m[2][3]/m[0][0] << endl;
        }

    } catch (exception& e) {

        cerr << e.what() << endl;

    }

    return 0;
}
