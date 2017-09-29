//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: crossgrating.cpp
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
#include "crossgrating.h"
#include "crossgrating2.h"
#include "gcross.h"
#include "matrixmath.h"
#include "scatmech.h"
#include "fft.h"
#include "bobvlieg.h"
#include <algorithm>
 
using namespace std;

namespace SCATMECH {

    namespace {
        void print(CFARRAY a, int n1, int n2, ostream& os,int option=0)
        {
            a.array(n1,n2);
            for (int i=1; i<=n1; ++i) {
                for (int j=1; j<=n2; ++j) {
                    if (option==0) os << abs(a(i,j));
                    else if (option==1) os << real(a(i,j));
                    else if (option==2) os << imag(a(i,j));
                    if (j!=n2) os << tab;
                }
                os << endl;
            }
        }
    }

    void CrossGrating::setup()
    {
        Model::setup();
    }

    void Gridded_CrossGrating::FourierFactorize()
    {
		using BobVlieg_Supp::mpow;
        CrossGrating::zeta = Gridded_CrossGrating::zeta;
        CrossGrating::d1 = Gridded_CrossGrating::d1;
        CrossGrating::d2 = Gridded_CrossGrating::d2;

        int i,j,k,l;
        int level;

        int M1 = 2*order1+1;
        int M2 = 2*order2+1;

        int M1a = 4*order1+1;
        int M2a = 4*order2+1;

        int MM = M1*M2;
        int N1 = grid1;
        int N2 = grid2;

        E0.allocate(M1,M2,M1,M2,levels);
        E1.allocate(M1,M2,M1,M2,levels);
        E2.allocate(M1,M2,M1,M2,levels);
        E3.allocate(M1,M2,M1,M2,levels);

        CFARRAY epsFT0(M1a,M2a);
        CFARRAY epsFT1(M1a,M2a);
        CFARRAY epsFT2(M1,M1,M2a);
        CFARRAY epsFT3(M2,M2,M1a);

        CFARRAY epsFT0a(M1a,N2);
        CFARRAY epsFT1a(M1a,N2);
        CFARRAY epsFT2a(M1a,N2);
        CFARRAY epsFT3a(M2a,N1);

        CFARRAY epsFT2b(M1,M1,N2);
        CFARRAY epsFT3b(M2,M2,N1);

        CFARRAY matrix1(M1,M1);
        CFARRAY matrix2(M2,M2);
        CFARRAY matrix3(M1,M2,M1,M2);

        CFARRAY tempN1(N1,1);
        CFARRAY tempN2(N2,1);

        //
        // Create FTs ...
        //
        // For each level in the stack...
        for (level=1; level<=levels; ++level) {

            // Start with FT'ing along 1st dimension...
            for (j=1; j<=N2; ++j) {
                for (i=1; i<=N1; ++i) tempN1(i) = eps(i,j,level);
                fft1d(tempN1,N1,-1);
                for (i=1; i<=M1a; ++i) epsFT0a(i,j) = tempN1((i-M1+N1)%N1+1)/(double)N1;

                for (i=1; i<=N1; ++i) tempN1(i) = 1./eps(i,j,level);
                fft1d(tempN1,N1,-1);
                for (i=1; i<=M1a; ++i) {
                    epsFT1a(i,j) = tempN1((i-M1+N1)%N1+1)/(double)N1;
                    epsFT2a(i,j) = tempN1((i-M1+N1)%N1+1)/(double)N1;
                }
            }

            // First FT for type 3 is along x2 direction...
            for (j=1; j<=N1; ++j) {
                for (i=1; i<=N2; ++i) tempN2(i) = 1./eps(j,i,level);
                fft1d(tempN2,N2,-1);
                for (i=1; i<=M2a; ++i) epsFT3a(i,j) = tempN2((i-M2+N2)%N2+1)/(double)N2;
            }

            // For types 0 and 1, finish the FT along the other dimension...
            for (i=1; i<=M1a; ++i) {
                for (j=1; j<=N2; ++j) tempN2(j) = epsFT0a(i,j);
                fft1d(tempN2,N2,-1);
                for (j=1; j<=M2a; ++j) epsFT0(i,j) = tempN2((j-M2+N2)%N2+1)/(double)N2;

                for (j=1; j<=N2; ++j) tempN2(j) = epsFT1a(i,j);
                fft1d(tempN2,N2,-1);
                for (j=1; j<=M2a; ++j) epsFT1(i,j) = tempN2((j-M2+N2)%N2+1)/(double)N2;
            }

            // Invert the partial 1/e matrices for type 2 ...
            for (i=1; i<=N2; ++i) {
                for (j=1; j<=M1; ++j) {
                    for (k=1; k<=M1; ++k) {
                        matrix1(j,k) = epsFT2a(j-k+2*order1+1,i);
                    }
                }
                Inverse(matrix1,M1);
                for (j=1; j<=M1; ++j) {
                    for (k=1; k<=M1; ++k) {
                        epsFT2b(j,k,i) = matrix1(j,k);
                    }
                }
            }

            // Invert the partial 1/e matrices for type 3 ...
            for (i=1; i<=N1; ++i) {
                for (j=1; j<=M2; ++j) {
                    for (k=1; k<=M2; ++k) {
                        matrix2(j,k) = epsFT3a(j-k+2*order2+1,i);
                    }
                }
                Inverse(matrix2,M2);
                for (j=1; j<=M2; ++j) {
                    for (k=1; k<=M2; ++k) {
                        epsFT3b(j,k,i) = matrix2(j,k);
                    }
                }
            }

            // Then finish FT'ing the type 2 matrices by integrating along x2 direction...
            for (i=1; i<=M1; ++i) {
                for (j=1; j<=M1; ++j) {
                    for (l=1; l<=N2; ++l) tempN2(l) = epsFT2b(i,j,l);
                    fft1d(tempN2,N2,-1);
                    for (l=1; l<=M2a; ++l) epsFT2(i,j,l) = tempN2((l-M2+N2)%N2+1)/(double)N2;
                }
            }

            // Then finish FT'ing the type 3 matrices by integrating along x1 direction...
            for (i=1; i<=M2; ++i) {
                for (j=1; j<=M2; ++j) {
                    for (l=1; l<=N1; ++l) tempN1(l) = epsFT3b(i,j,l);
                    fft1d(tempN1,N1,-1);
                    for (l=1; l<=M1a; ++l) epsFT3(i,j,l) = tempN1((l-M1+N1)%N1+1)/(double)N1;
                }
            }

            for (i=1; i<=M1; ++i) {
                for (j=1; j<=M2; ++j) {
                    for (k=1; k<=M1; ++k) {
                        for (l=1; l<=M2; ++l) {
                            E2(i,j,k,l,level) = epsFT2(i,k,j-l+2*order2+1);
                            E3(i,j,k,l,level) = epsFT3(j,l,i-k+2*order1+1);
                        }
                    }
                }
            }

            for (i=1; i<=M1; ++i) {
                for (j=1; j<=M2; ++j) {
                    for (k=1; k<=M1; ++k) {
                        for (l=1; l<=M2; ++l) {
                            matrix3(i,j,k,l) = epsFT0(i-k+2*order1+1,j-l+2*order2+1);
                        }
                    }
                }
            }

            Inverse(matrix3,MM);

            for (i=1; i<=M1; ++i) {
                for (j=1; j<=M2; ++j) {
                    for (k=1; k<=M1; ++k) {
                        for (l=1; l<=M2; ++l) {
                            E0(i,j,k,l,level) = matrix3(i,j,k,l);
                        }
                    }
                }
            }

            for (i=1; i<=M1; ++i) {
                for (j=1; j<=M2; ++j) {
                    for (k=1; k<=M1; ++k) {
                        for (l=1; l<=M2; ++l) {
                            matrix3(i,j,k,l) = epsFT1(i-k+2*order1+1,j-l+2*order2+1);
                        }
                    }
                }
            }

            Inverse(matrix3,MM);

            for (i=1; i<=M1; ++i) {
                for (j=1; j<=M2; ++j) {
                    for (k=1; k<=M1; ++k) {
                        for (l=1; l<=M2; ++l) {
                            E1(i,j,k,l,level) = matrix3(i,j,k,l);
                        }
                    }
                }
            }
        }
    }

    void CrossGrating::set_parameter_base(
        const STRING& parameter, ///< The parameter name
        const STRING& value      ///< String represention of a value
    )
    {
        if (parameter=="lambda") {
            lambda = from_string<double>(value);
            set_recalc(0xFF);
        } else if (parameter=="order1") {
            order1 = from_string<int>(value);
            set_recalc(0xFF);
        } else if (parameter=="order2") {
            order2 = from_string<int>(value);
            set_recalc(0xFF);
        } else if (parameter=="SaveGratingE0") {
            SETUP();
            int M1 = 2*order1+1;
            int M2 = 2*order2+1;

            ofstream image(value.c_str());
            if (!image) error("Cannot open file: " + value);
            for (int level=1; level<=levels; ++level) {
                print(E0(1,1,1,1,level),M1*M2,M1*M2,image);
                image << endl;
            }
            set_recalc(0);
        } else if (parameter=="SaveGratingE1") {
            SETUP();
            int M1 = 2*order1+1;
            int M2 = 2*order2+1;

            ofstream image(value.c_str());
            if (!image) error("Cannot open file: " + value);
            for (int level=1; level<=levels; ++level) {
                print(E1(1,1,1,1,level),M1*M2,M1*M2,image);
                image << endl;
            }
            set_recalc(0);
        } else if (parameter=="SaveGratingE2") {
            SETUP();
            int M1 = 2*order1+1;
            int M2 = 2*order2+1;

            ofstream image(value.c_str());
            if (!image) error("Cannot open file: " + value);
            for (int level=1; level<=levels; ++level) {
                print(E2(1,1,1,1,level),M1*M2,M1*M2,image);
                image << endl;
            }
            set_recalc(0);
        } else if (parameter=="SaveGratingE3") {
            SETUP();
            int M1 = 2*order1+1;
            int M2 = 2*order2+1;

            ofstream image(value.c_str());
            if (!image) error("Cannot open file: " + value);
            for (int level=1; level<=levels; ++level) {
                print(E3(1,1,1,1,level),M1*M2,M1*M2,image);
                image << endl;
            }
            set_recalc(0);
        } else if (parameter=="SaveLevelThickness") {
            SETUP();
            ofstream file(value.c_str());
            if (!file) error("Cannot open file: " + value);
            for (int level=1; level<=levels; ++level) {
                file << level << tab << thick(level) << endl;
            }
        } else {
            Model::set_parameter_base(parameter,value);
        }
    }

    void Gridded_CrossGrating::set_parameter_base(
        const STRING& parameter, ///< The parameter name
        const STRING& value      ///< String represention of a value
    )
    {
        if (parameter=="SaveGratingImage") {
            SETUP();
            ofstream image(value.c_str());
            if (!image) error("Cannot open file: " + value);
            for (int level=1; level<=levels; ++level) {
                print(eps(1,1,level),grid1,grid2,image);
                image << endl;
            }
            set_recalc(0);
        } else if (parameter=="SaveGratingImageReal") {
            SETUP();
            ofstream image(value.c_str());
            if (!image) error("Cannot open file: " + value);
            for (int level=1; level<=levels; ++level) {
                print(eps(1,1,level),grid1,grid2,image,1);
                image << endl;
            }
            set_recalc(0);
        } else if (parameter=="SaveGratingImageImag") {
            SETUP();
            ofstream image(value.c_str());
            if (!image) error("Cannot open file: " + value);
            for (int level=1; level<=levels; ++level) {
                print(eps(1,1,level),grid1,grid2,image,2);
                image << endl;
            }
            set_recalc(0);
        } else {
            CrossGrating::set_parameter_base(parameter,value);
        }
    }

    void Gridded_CrossGrating::getxy(int i,int j,double &x, double &y)
    {
        double x1 = (i-0.5-grid1/2.)/grid1*d1;
        double x2 = (j-0.5-grid2/2.)/grid2*d2;

        x = x1 + x2*sinzeta;
        y = x2*coszeta;
    }

    void Register(const CrossGrating* x)
    {
        static bool regd=false;
        if (!regd) {
            regd=true;

            Register_Model(CrossGrating);
            Register_Model(OneD_CrossGrating);
            Register_Model(Overlaid_CrossGrating);
            Register_Model(Overlaid_1D_CrossGrating);
            Register_Model(Null_CrossGrating);
            Register_Model(Gridded_CrossGrating);
            Register_Model(Cylinder_CrossGrating);
            Register_Model(Sphere_CrossGrating);
            Register_Model(Pyramidal_Pit_CrossGrating);
            Register_Model(Generic_CrossGrating);
        }
    }

    DEFINE_VIRTUAL_MODEL(CrossGrating,Model,"Generalized cross grating");
    DEFINE_PARAMETER(CrossGrating,dielectric_function,medium_i,"Incident medium","(1,0)",0xFF);
    DEFINE_PARAMETER(CrossGrating,dielectric_function,medium_t,"Transmitted medium","(4.05,0.05)",0xFF);

    DEFINE_VIRTUAL_MODEL(Gridded_CrossGrating,CrossGrating,"Cross grating defined on a discrete grid");
    DEFINE_PARAMETER(Gridded_CrossGrating,double,zeta,"Angle of lattice vectors from perpendicular [deg]","0.",0xFF);
    DEFINE_PARAMETER(Gridded_CrossGrating,double,d1,"Lattice constant #1 [um]","0.5",0xFF);
    DEFINE_PARAMETER(Gridded_CrossGrating,double,d2,"Lattice constant #2 [um]","0.5",0xFF);
    DEFINE_PARAMETER(Gridded_CrossGrating,int,grid1,"Number of sampling points in direction #1","1024",0xFF);
    DEFINE_PARAMETER(Gridded_CrossGrating,int,grid2,"Number of sampling points in direction #2","1024",0xFF);

}
