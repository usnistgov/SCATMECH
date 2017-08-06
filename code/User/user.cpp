//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: user.cpp
//**
//** Thomas A. Germer
//** Optical Technology Division,
//** National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//** Version: 7.00 (January 2015)
//**
//******************************************************************************
//
// TODO: Rename user.cpp and user.h to something more appropriate, and change the
//       following include file to the respective include file.
//
#include "user.h"

using namespace std;
using namespace SCATMECH;

//
// TODO: Search and replace "User_BRDF_Model" with an appropriate name.
// TODO: Search and replace "Parent_BRDF_Model" with the appropriate parent model.
//

void
User_BRDF_Model::
setup()
{
    Parent_BRDF_Model::setup();
    // TODO: Add any necessary operations here that need to be performed
    //       whenever parameters change...
    // ...

}

//
// TODO: Decide whether this model is a JonesMatrix model or a MuellerMatrix
//       model and de-comment one of the following two functions:
//       If the model inherits Local_BRDF_Model, then these functions should
//       be renamed jonesDSC or muellerDSC, respectively.
//
//	 The variables thetai, thetas, phis, and rotation are available to
//       these routines, and designate the incident and viewing directions in radians
/*
JonesMatrix
User_BRDF_Model::
jones()
{
    // This macro checks the variable recalc and calls setup if necc.
    SETUP(); 

    JonesMatrix result;
    //
    // TODO: Code for model goes here...
    //

    return result;
}
*/

/*
MuellerMatrix
User_BRDF_Model::
mueller()
{
    // This macro checks the variable recalc and calls setup if necc.
    SETUP(); 

    MuellerMatrix result;
    //
    // TODO: Code for model goes here...
    //

    return result;
}
*/

//
// TODO: Give the model a description
//
DEFINE_MODEL(User_BRDF_Model,
	     Parent_BRDF_Model,
	     "Description");

//
// TODO: Use the DEFINE_PARAMETER macro for each model parameter.
//       See inherit.h for details...
//
DEFINE_PARAMETER(User_BRDF_Model, 
		 double,         // The data type
		 shininess,      // The parameter name
		 "The material shininess", // A description
		 "default");     // A string version of the default value

//
// TODO: A call to Register_Model(User_BRDF_Model) must be made before
//       User_BRDF_Model is known to library.
//       Add this function call (it's actually a macro) to either main()
//       or some function called by main().
//

