//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: user.h
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
//TODO: Change the following macro names to something appropriate...
//
#ifndef SCATMECH_USER_H
#define SCATMECH_USER_H

//
// TODO: Search and replace "User_BRDF_Model" with an appropriate name.
// TODO: Search and replace "Parent_BRDF_Model" with the appropriate parent model.
//       You may need to use the qualifier SCATMECH:: on the parent class
//

#include "brdf.h"

class User_BRDF_Model : public Parent_BRDF_Model
{
protected:
    // The function setup() allows one to perform one-time operations whenever
    // member data changes:
    virtual void setup();

    //
    // TODO: Decide whether this model is a JonesMatrix model or a MuellerMatrix
    //       model and de-comment one of the following two functions.
    //       If the model is a Local_BRDF_Model, change jones to jonesDSC and
    //       mueller to muellerDSC.
    //
    // virtual SCATMECH::JonesMatrix jones();
    // virtual SCATMECH::MuellerMatrix mueller();

private:

    DECLARE_MODEL();

    //
    // TODO: Replace the following parameter declaration with declarations
    //       appropriate for the model...
    //
    DECLARE_PARAMETER(double,shininess);
};

#endif

