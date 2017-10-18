//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: gcross.h
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//******************************************************************************
#ifndef SCATMECH_GCROSS_H
#define SCATMECH_GCROSS_H

#include "vector3d.h"
#include "matrix3d.h"
#include "crossgrating.h"
#include "crossgrating2.h"

namespace SCATMECH {

    class Generic_CrossGrating : public Gridded_CrossGrating {
        public:

            DECLARE_MODEL();
            DECLARE_PARAMETER(std::string,filename);
            DECLARE_PARAMETER(std::string,pstring);
            DECLARE_PARAMETER(int,nlayers);

        public:

            struct boundary {
                int n;
                std::string *text;
                Vector vertex1,vertex2,vertex3;
                COMPLEX mat1,mat2;
            };
            typedef std::vector<boundary> boundary_vector;

            struct bounds {
                double x;
                COMPLEX mat1,mat2;
                int n; // The line number of this boundary
                std::string *text;

                bool operator<(const bounds& b) {
                    if (x==b.x) return mat2==b.mat1;
                    else return x<b.x;
                }
                bool operator==(const bounds& b) {
                    if (x!=b.x) return false;
                    if (mat1!=b.mat1) return false;
                    if (mat2!=b.mat2) return false;
                    return true;
                }

            };

            const boundary_vector& Get_Boundaries() {
                SETUP();
                return boundaries;
            }

            void print_parameters(std::ostream& os,const std::string& prefix) const;

        protected:
            void setup();

            // The following three override the same named function of Model...
            virtual STRING get_parameter_base(const STRING& parameter) const;
            virtual void set_parameter_base(const STRING& parameter, const STRING& value);

        private:
            std::vector<std::vector<std::string> > material;
            std::vector<std::vector<double> > position;
            std::vector<double> thickness;

            std::string filecontents;
            std::string old_filename;

            typedef std::map<std::string,COMPLEX> epsmap;
            epsmap eps0;

            typedef std::map<std::string,double> varsmap;
            varsmap vars;

            typedef std::map<std::string,Vector> verticesmap;
            verticesmap vertices;


            COMPLEX eps_t;
            COMPLEX eps_i;

            varsmap parameters;
            varsmap override;

            boundary_vector boundaries;

            COMPLEX find_epsilon(const Vector& s,bool throwerror);
            Vector get_triplet_value(std::istream& is) const;
            void error(const std::string& message) const;
    };

}

#endif
