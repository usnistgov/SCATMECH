//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: firstdiffuse.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_FIRSTDIFFUSE_H
#define SCATMECH_FIRSTDIFFUSE_H

#include "brdf.h"
#include "askuser.h"
#include "filmtran.h"


namespace SCATMECH {

    class Phase_Function : public Model
    {
        public:
            Phase_Function() {}
            virtual double f(double theta)=0;
            DECLARE_MODEL();
    };

    typedef Model_Ptr<Phase_Function> Phase_Function_Ptr;

    class First_Diffuse_BRDF_Model : public BRDF_Model
    {
        public:
            First_Diffuse_BRDF_Model();

            DECLARE_MODEL();
            DECLARE_PARAMETER(double,l_scat);
            DECLARE_PARAMETER(Phase_Function_Ptr,phase_function);
            DECLARE_PARAMETER(StackModel_Ptr,stack);

        protected:
            virtual void setup();
            virtual JonesMatrix jones();
    };


    class Henyey_Greenstein_Phase_Function: public Phase_Function
    {
        public:
            virtual double f(double theta);
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,g);
    };

    class Double_Henyey_Greenstein_Phase_Function: public Phase_Function
    {
        public:
            virtual double f(double theta);
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,g);
            DECLARE_PARAMETER(double,c);
    };

    class Isotropic_Phase_Function: public Phase_Function
    {
        public:
            virtual double f(double theta) {
                return 0.25/pi;
            }
            DECLARE_MODEL();
    };

    class Rayleigh_Phase_Function: public Phase_Function
    {
        public:
            virtual double f(double theta) {
                return 0.75*(1.+sqr(cos(theta)));
            }
            DECLARE_MODEL();
    };

    class Table_Phase_Function: public Phase_Function
    {
        public:
            virtual double f(double theta);
            DECLARE_MODEL();
            DECLARE_PARAMETER(Table,table);
    };

    class Reynolds_McCormick_Phase_Function: public Phase_Function
    {
        public:
            virtual double f(double theta);
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,g);
            DECLARE_PARAMETER(double,alpha);
        protected:
            void setup();
        private:
            double K;
    };

    class Double_Reynolds_McCormick_Phase_Function: public Phase_Function
    {
        public:
            virtual double f(double theta);
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,g);
            DECLARE_PARAMETER(double,alpha);
            DECLARE_PARAMETER(double,c);
        protected:
            void setup();
        private:
            double K;
    };

    class Legendre_Phase_Function: public Phase_Function
    {
        public:
            virtual double f(double theta);
            DECLARE_MODEL();
            DECLARE_PARAMETER(double,c0);
            DECLARE_PARAMETER(double,c1);
            DECLARE_PARAMETER(double,c2);
            DECLARE_PARAMETER(double,c3);
            DECLARE_PARAMETER(double,c4);
            DECLARE_PARAMETER(double,c5);
    };


    void Register(const First_Diffuse_BRDF_Model*);
    void Register(const Phase_Function* x);


}

#endif

