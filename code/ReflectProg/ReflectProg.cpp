//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: ReflectProg.cpp
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
#include "filmtran.h"

using namespace std;    // Use unqualified names for Standard C++ library
using namespace SCATMECH;


int main()
{
    try {

        // Query user for wavelength...
        double lambda = AskUser("Wavelength",0.633);

        // Query user for substrate...
        dielectric_function substrate =
            dielectric_function::AskUser("Substrate",optical_constant(1.5,0));

        // Query user for films...
        dielectric_stack stack = dielectric_stack::AskUser("Stack","");

        // Loop through angles...
        for (double theta=0; theta<90; theta+=1) {

            // Print out s- and p-polarized reflectances...
            cout << theta << tab
                 << norm(stack.rp12(theta*deg,lambda,vacuum,substrate)) << tab
                 << norm(stack.rs12(theta*deg,lambda,vacuum,substrate)) << endl;
        }

    } catch (exception& e) {

        cerr << e.what() << endl;

    }

    return 0;
}
