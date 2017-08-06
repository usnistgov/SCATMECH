//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: RCWProg.cpp
//**
//** Thomas A. Germer
//** Optical Technology Division,
//** National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//** Version 7.00 (January 2015)
//**
//******************************************************************************
#include "rcw.h"

using namespace std;
using namespace SCATMECH;

int main()
{
    try {

        // Create object...
        RCW_Model RCW;

        // Query user for parameters...
        RCW.AskUser();

        // This keeps RCW from being verbose...
        RCW.set_quiet(1);

        // Display all of the available parameters...
        cerr << endl << "Available parameters:" << endl << endl;
        RCW.print_parameters(cerr);
        cerr << endl;

        // Get parameter to vary, its limits, and its range...
        string parameter=AskUser("Parameter to vary",string("lambda"));
        double begin=AskUser("Initial value",0.200);
        double end=AskUser("Final value",0.800);
        double step=AskUser("Step value",0.005);

        // Output file header...
        cout << parameter << tab << "alpha" << tab << "beta" << endl;

        // Scan parameter...
        for (double value = begin; value <= end; value += step) {

            // Set the parameter...
            RCW.set_parameter(parameter,value);

            // Get the Mueller matrix intensity of the specular (0 order) reflection...
            MuellerMatrix M = RCW.GetIntensity(0);

            // Get the Stokes vector for 45 deg polarized incident light...
            StokesVector S = M*StokesVectorUnitLinear(45.*deg);

            // Output the normalized Stokes parameters measured by a rotating analyzer ellipsometer...
            cout << value << tab << S[1]/S[0] << tab << S[2]/S[0] << endl;
        }

    } catch (exception& e) {

        cerr << e.what() << endl;

    }
    return 0;
}
