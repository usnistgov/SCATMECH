//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: allrough.h
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
#ifndef SCATMECH_ALLROUGH_H
#define SCATMECH_ALLROUGH_H

#include "roughnes.h"


namespace SCATMECH {


    class Correlated_Roughness_Stack_BRDF_Model : public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(dielectric_stack,stack);

        protected:

            void setup();
            JonesMatrix jones();

        private:

            Roughness_Stack_BRDF_Model model;
    };

    class Uncorrelated_Roughness_Stack_BRDF_Model : public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(dielectric_stack,stack);

        protected:

            void setup();
            MuellerMatrix mueller();

        private:

            Roughness_Stack_BRDF_Model model;
    };

    class Growth_Roughness_Stack_BRDF_Model : public Roughness_BRDF_Model
    {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(dielectric_stack,stack);
            DECLARE_PARAMETER(double,relaxation);
            DECLARE_PARAMETER(double,exponent);
            DECLARE_PARAMETER(PSD_Function_Ptr,intrinsic);

        protected:

            void setup();
            MuellerMatrix mueller();

        private:

            Roughness_Stack_BRDF_Model model;
    };

} // namespace SCATMECH;


#endif
