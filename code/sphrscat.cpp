//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: sphrscat.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#include "sphrscat.h"
#include "askuser.h"
#include "matrix3d.h"

using namespace std;


namespace SCATMECH {

    DEFINE_VIRTUAL_MODEL(Free_Space_Scatterer,Model,
                         "Generalized free-space scatterer");

    DEFINE_PARAMETER(Free_Space_Scatterer,double,lambda,"Wavelength [um]","0.532",0xFF);
    DEFINE_PARAMETER(Free_Space_Scatterer,dielectric_function,medium,"Optical properties of the surrounding medium","(1,0)",0xFF);

    DEFINE_VIRTUAL_MODEL(SphericalScatterer,Free_Space_Scatterer,
                         "A spherically-symmetric free-space scatterer");

    DEFINE_PARAMETER(SphericalScatterer,double,radius,"Radius [um]","0.05",0xFF);
    DEFINE_PARAMETER(SphericalScatterer,dielectric_function,sphere,"Optical properties of the sphere","(1.59,0)",0xFF);


} // namespace SCATMECH



