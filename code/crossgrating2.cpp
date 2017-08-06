//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: crossgrating2.cpp
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
#include "crossgrating2.h"

using namespace std;


namespace SCATMECH {

    void Cylinder_CrossGrating::setup()
    {
        Gridded_CrossGrating::setup();

        levels = nlevels;
        eps.allocate(grid1,grid2,levels);
        thick.allocate(levels);

        COMPLEX einside = epsilon(inside);
        COMPLEX eoutside = epsilon(outside);

        for (int n=1; n<=levels; ++n) {
            thick(n) = thickness/levels;
        }
        for (int i=1; i<=grid1; ++i) {
            for (int j=1; j<=grid2; ++j) {
                double x,y;
                getxy(i,j,x,y);

                double theta=atan2(y,x);
                double r = sqrt(x*x+y*y);

                double _rtop = rtop.value(theta/deg);
                double _rbottom = rbottom.value(theta/deg);

                for (int n=1; n<=levels; ++n) {
                    double radius = (n-0.5)/levels*(_rtop-_rbottom)+_rbottom;

                    if (r<=radius) {
                        eps(i,j,n) = einside;
                    } else {
                        eps(i,j,n) = eoutside;
                    }
                }
            }
        }

        FourierFactorize();
    }

    void Sphere_CrossGrating::setup()
    {
        Gridded_CrossGrating::setup();

        if (above<0) error("abovelevels<0");
        if (below<0) error("belowlevels<0");

        int abovelevels = (above!=0. ? 1 : 0);
        int belowlevels = (below!=0. ? 1 : 0);
        levels =  nlevels + abovelevels + belowlevels ;

        eps.allocate(grid1,grid2,levels);
        thick.allocate(levels);

        COMPLEX esurround = (COMPLEX)surrounding.epsilon(lambda);
        COMPLEX esphere = (COMPLEX)sphere.epsilon(lambda);
        for (int i=1; i<=grid1*grid2*levels; ++i) eps(i)=esurround;

        if (below!=0.) thick(1) = below;
        if (above!=0.) thick(levels) = above;


        double radius = diameter/2.;
        double t = diameter/nlevels;

        for (int n=1; n<=nlevels; ++n) {

            double r = (n<=(nlevels+1)/2) ? (2*n-1)*diameter/nlevels/2 : (2*(nlevels-n+1)-1)*diameter/nlevels/2;

            if (n==1 || n==nlevels) {
                thick(n+belowlevels)=radius-sqrt(sqr(radius)-sqr(r+radius/nlevels));
            } else if ((double)n==(nlevels+1.)/2.) {
                thick(n+belowlevels)=2*sqrt(sqr(radius)-sqr(r-radius/nlevels));
            } else {
                thick(n+belowlevels)=sqrt(sqr(radius)-sqr(r-radius/nlevels))-sqrt(sqr(radius)-sqr(r+radius/nlevels));
            }

            for (int i=1; i<=grid1; ++i) {
                for (int j=1; j<=grid2; ++j) {
                    double x,y;
                    getxy(i,j,x,y);

                    if (sqr(x)+sqr(y)<sqr(r)) {
                        eps(i,j,n+belowlevels) = esphere;
                    }
                }
            }
        }

        FourierFactorize();
    }

    void Pyramidal_Pit_CrossGrating::setup()
    {
        Gridded_CrossGrating::setup();

        levels = nlevels;

        eps.allocate(grid1,grid2,nlevels);
        thick.allocate(levels);

        COMPLEX medium_t_eps = (COMPLEX)medium_t.epsilon(lambda);
        COMPLEX medium_i_eps = (COMPLEX)medium_i.epsilon(lambda);
        for (int i=1; i<=grid1*grid2*levels; ++i) eps(i)=medium_t_eps;

        for (int n=1; n<=levels; ++n) {
            thick(n) = depth/levels;
            double hside = side*(n-0.5)/levels/2.;
            for (int i=1; i<=grid1; ++i) {
                for (int j=1; j<=grid2; ++j) {
                    double x,y;
                    getxy(i,j,x,y);

                    if (fabs(x)<hside && fabs(y)<hside) eps(i,j,n) = medium_i_eps;
                }
            }
        }
        FourierFactorize();
    }

    //
    // The routine FillE0E1E2E3Matrices fills the E-matrices using 1D Grating data.
    //
    // "level2d" is the level number of the 2D grating.
    // "level1d" is the level number of the 1D grating.
    // "flip" specifies
    // (if false) the grating direction is aligned along the "1" axis, or
    // (if true) the grating direction is aligned along the "2" axis.
    //
    //
    void FillE0E1E2E3Matrices(CFARRAY& E0,CFARRAY& E1,CFARRAY& E2,CFARRAY& E3, int M1,int M2, int level2d, int level1d, Grating_Ptr& grating, bool flip)
    {
        int M = flip ? M2 : M1;

        CFARRAY matrix(M,M);
        for (int i=1; i<=M; ++i) {
            for (int k=1; k<=M; ++k) {
                matrix(i,k) = (grating->fourierz(i-k,level1d,0));
            }
        }
        Inverse(matrix,M);
        for (int i=1; i<=M1; ++i) {
            for (int j=1; j<=M2; ++j) {
                for (int k=1; k<=M1; ++k) {
                    for (int l=1; l<=M2; ++l) {
                        if (flip) {
                            E0(i,j,k,l,level2d) = (i==k) ? matrix(j,l) : 0.;
                        } else {
                            E0(i,j,k,l,level2d) = (j==l) ? matrix(i,k) : 0.;
                        }
                    }
                }
            }
        }

        for (int i=1; i<=M; ++i) {
            for (int j=1; j<=M; ++j) {
                matrix(i,j) = (grating->fourierx(i-j,level1d,1));
            }
        }
        Inverse(matrix,M);
        for (int i=1; i<=M1; ++i) {
            for (int j=1; j<=M2; ++j) {
                for (int k=1; k<=M1; ++k) {
                    for (int l=1; l<=M2; ++l) {
                        if (flip) {
                            E1(i,j,k,l,level2d) = (i==k) ? matrix(j,l) : 0.;
                            E3(i,j,k,l,level2d) = (i==k) ? matrix(j,l) : 0.;
                        } else {
                            E1(i,j,k,l,level2d) = (j==l) ? matrix(i,k) : 0.;
                            E2(i,j,k,l,level2d) = (j==l) ? matrix(i,k) : 0.;
                        }
                    }
                }
            }
        }

        for (int i=1; i<=M; ++i) {
            for (int k=1; k<=M; ++k) {
                matrix(i,k) = (grating->fouriery(i-k,level1d,0));
            }
        }
        for (int i=1; i<=M1; ++i) {
            for (int j=1; j<=M2; ++j) {
                for (int k=1; k<=M1; ++k) {
                    for (int l=1; l<=M2; ++l) {
                        if (flip) {
                            E2(i,j,k,l,level2d) = (i==k) ? matrix(j,l) : 0.;
                        } else {
                            E3(i,j,k,l,level2d) = (j==l) ? matrix(i,k) : 0.;
                        }
                    }
                }
            }
        }
    }


    void OneD_CrossGrating::setup()
    {
        CrossGrating::setup();

        if (grating->get_lambda()!=lambda) grating->set_lambda(lambda);
        if ((COMPLEX)grating->get_medium_i().index(lambda)!=(COMPLEX)medium_i.index(lambda)) error("grating.medium_i!=medium_i");
        if ((COMPLEX)grating->get_medium_t().index(lambda)!=(COMPLEX)medium_t.index(lambda)) error("grating.medium_t!=medium_t");
        if (grating->is_magnetic()) error("Magnetic materials not supported");

        levels = grating->get_levels();

        double period = grating->get_period();
        d1 = period/cos(zeta*deg);
        CrossGrating::d2 = d2;
        CrossGrating::zeta = zeta;

        thick.allocate(levels);

        int M1 = 2*order1+1;
        int M2 = 2*order2+1;

        E0.allocate(M1,M2,M1,M2,levels);
        E1.allocate(M1,M2,M1,M2,levels);
        E2.allocate(M1,M2,M1,M2,levels);
        E3.allocate(M1,M2,M1,M2,levels);

        for (int level=1; level<=levels; ++level) {
            thick(level) = grating->get_thickness(levels-level);
            FillE0E1E2E3Matrices(E0,E1,E2,E3,M1,M2,level,levels-level,grating,false);
        }
    }

    void Overlaid_CrossGrating::setup()
    {
        CrossGrating::setup();

        top->set_lambda(lambda);
        top->set_order1(order1);
        top->set_order2(order2);

        bottom->set_lambda(lambda);
        bottom->set_order1(order1);
        bottom->set_order2(order2);

        if (top->get_zeta()!=bottom->get_zeta()) error("top.zeta != bottom.zeta");
        if (top->get_d1()!=bottom->get_d1()) error("top.d1 != bottom.d1");
        if (top->get_d2()!=bottom->get_d2()) error("top.d2 != bottom.d2");
        if (top->get_medium_i().index(lambda)!=medium_i.index(lambda)) error("top.medium_i!=medium_i");
        if (bottom->get_medium_t().index(lambda)!=medium_t.index(lambda)) error("bottom.medium_t!=medium_t");
        if (top->get_medium_t().index(lambda)!=bottom->get_medium_i().index(lambda)) error("top.medium_t!=bottom.medium_i");

        int tlevels = top->get_levels();
        int blevels = bottom->get_levels();

        levels = tlevels + blevels;
        if (separation>0.) ++levels;

        d1 = top->get_d1();
        d2 = top->get_d2();
        zeta = top->get_zeta();

        CFARRAY tE0 = top->get_E0();
        CFARRAY tE1 = top->get_E1();
        CFARRAY tE2 = top->get_E2();
        CFARRAY tE3 = top->get_E3();

        CFARRAY bE0 = bottom->get_E0();
        CFARRAY bE1 = bottom->get_E1();
        CFARRAY bE2 = bottom->get_E2();
        CFARRAY bE3 = bottom->get_E3();

        int M1 = 2*order1+1;
        int M2 = 2*order2+1;

        E0.allocate(M1,M2,M1,M2,levels);
        E1.allocate(M1,M2,M1,M2,levels);
        E2.allocate(M1,M2,M1,M2,levels);
        E3.allocate(M1,M2,M1,M2,levels);
        thick.allocate(levels);

        double k1 = 2*pi/d1;
        double k2 = 2*pi/d2;
        COMPLEX cI(0,1);

        for (int level=1; level<=levels; ++level) {
            if (level<=blevels) {
                int blevel = level;
                thick(level) = bottom->get_thick(blevel);
                for (int i=1; i<=M1; ++i) {
                    for (int j=1; j<=M2; ++j) {
                        for (int k=1; k<=M1; ++k) {
                            for (int l=1; l<=M2; ++l) {
                                E0(i,j,k,l,level) = bE0(i,j,k,l,blevel);
                                E1(i,j,k,l,level) = bE1(i,j,k,l,blevel);
                                E2(i,j,k,l,level) = bE2(i,j,k,l,blevel);
                                E3(i,j,k,l,level) = bE3(i,j,k,l,blevel);
                            }
                        }
                    }
                }
            }
            if (level==blevels+1 && separation>0.) {
                thick(level) = separation;
                COMPLEX smedium = top->get_medium_t().epsilon(lambda);
                for (int i=1; i<=M1; ++i) {
                    for (int j=1; j<=M2; ++j) {
                        for (int k=1; k<=M1; ++k) {
                            for (int l=1; l<=M2; ++l) {
                                E0(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                                E1(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                                E2(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                                E3(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                            }
                        }
                    }
                }
            }
            int iseparation = separation>0.? 1 : 0;
            if (level>blevels + iseparation) {
                int tlevel = level - blevels - iseparation;
                thick(level) = top->get_thick(tlevel);
                for (int i=1; i<=M1; ++i) {
                    for (int j=1; j<=M2; ++j) {
                        for (int k=1; k<=M1; ++k) {
                            for (int l=1; l<=M2; ++l) {
                                COMPLEX phase = exp(-cI*(k1*overlay1*(i-k)+k2*overlay2*(j-l)));
                                E0(i,j,k,l,level) = tE0(i,j,k,l,tlevel)*phase;
                                E1(i,j,k,l,level) = tE1(i,j,k,l,tlevel)*phase;
                                E2(i,j,k,l,level) = tE2(i,j,k,l,tlevel)*phase;
                                E3(i,j,k,l,level) = tE3(i,j,k,l,tlevel)*phase;
                            }
                        }
                    }
                }
            }
        }
    }


    void Overlaid_1D_CrossGrating::setup()
    {
        CrossGrating::setup();

        top->set_lambda(lambda);

        bottom->set_lambda(lambda);

        if (top->get_medium_t().index(lambda)!=bottom->get_medium_i().index(lambda)) error("top.medium_t!=bottom.medium_i");
        if (top->get_medium_i().index(lambda)!=medium_i.index(lambda)) error("top.medium_i!=medium_i");
        if (bottom->get_medium_t().index(lambda)!=medium_t.index(lambda)) error("bottom.medium_t!=medium_t");
        if (angle == 0.) error("angle=0 not supported");

        zeta = 90.-angle;
        d1 = top->get_period()/sin(angle*deg);
        d2 = bottom->get_period()/sin(angle*deg);

        int tlevels = top->get_levels();
        int blevels = bottom->get_levels();

        levels = tlevels + blevels;
        if (separation>0.) ++levels;

        thick.allocate(levels);

        int M1 = 2*order1+1;
        int M2 = 2*order2+1;

        E0.allocate(M1,M2,M1,M2,levels);
        E1.allocate(M1,M2,M1,M2,levels);
        E2.allocate(M1,M2,M1,M2,levels);
        E3.allocate(M1,M2,M1,M2,levels);

        int iseparation = separation>0.? 1 : 0;

        for (int level=1; level<=levels; ++level) {
            if (level<=blevels) {
                int blevel = blevels-level;
                thick(level) = bottom->get_thickness(blevel);
                FillE0E1E2E3Matrices(E0,E1,E2,E3,M1,M2,level,blevel,bottom,true);
            } else if (level==blevels+1 && separation>0.) {
                thick(level) = separation;
                COMPLEX smedium = top->get_medium_t().epsilon(lambda);
                for (int i=1; i<=M1; ++i) {
                    for (int j=1; j<=M2; ++j) {
                        for (int k=1; k<=M1; ++k) {
                            for (int l=1; l<=M2; ++l) {
                                E0(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                                E1(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                                E2(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                                E3(i,j,k,l,level) = (i==k && j==l) ? smedium : 0.;
                            }
                        }
                    }
                }
            } else if (level>blevels + iseparation) {
                int tlevel = blevels-level+blevels+iseparation;
                thick(level) = top->get_thickness(tlevel);
                FillE0E1E2E3Matrices(E0,E1,E2,E3,M1,M2,level,tlevel,top,false);
            }
        }
    }


    void Null_CrossGrating::setup()
    {
        CrossGrating::setup();

        levels = 0;

        CrossGrating::d1 = d1;
        CrossGrating::d2 = d2;
        CrossGrating::zeta = zeta;

        int M1 = 2*order1+1;
        int M2 = 2*order2+1;

        E0.allocate(M1,M2,M1,M2,0);
        E1.allocate(M1,M2,M1,M2,0);
        E2.allocate(M1,M2,M1,M2,0);
        E3.allocate(M1,M2,M1,M2,0);
        thick.allocate(levels);
    }



    DEFINE_MODEL(Cylinder_CrossGrating,Gridded_CrossGrating,"Contact holes in a 2-d array");
    DEFINE_PARAMETER(Cylinder_CrossGrating,Table,rtop,"Radius of top of holes [um] as function of angle [deg]","0.1",0xFF);
    DEFINE_PARAMETER(Cylinder_CrossGrating,Table,rbottom,"Radius of bottom of holes [um] as function of angle [deg]","0.1",0xFF);
    DEFINE_PARAMETER(Cylinder_CrossGrating,double,thickness,"Thickness of grating [um]","0.1",0xFF);
    DEFINE_PARAMETER(Cylinder_CrossGrating,int,nlevels,"Number of levels in grating","1",0xFF);
    DEFINE_PARAMETER(Cylinder_CrossGrating,dielectric_function,inside,"Medium inside holes","(1,0)",0xFF);
    DEFINE_PARAMETER(Cylinder_CrossGrating,dielectric_function,outside,"Medium outside holes","(1.5,0)",0xFF);

    DEFINE_MODEL(OneD_CrossGrating,CrossGrating,"One dimensional grating");
    DEFINE_PARAMETER(OneD_CrossGrating,double,d2,"Lattice constant #2 [um]","0.5",0xFF);
    DEFINE_PARAMETER(OneD_CrossGrating,double,zeta,"Angle of lattice vectors from perpendicular [deg]","0",0xFF);
    DEFINE_PTRPARAMETER(OneD_CrossGrating,Grating_Ptr,grating,"Grating","Single_Line_Grating",0xFF);

    DEFINE_MODEL(Overlaid_CrossGrating,CrossGrating,"A crossed grating on top of another");
    DEFINE_PTRPARAMETER(Overlaid_CrossGrating,CrossGrating_Ptr,top,"Top grating","OneD_CrossGrating",0xFF);
    DEFINE_PTRPARAMETER(Overlaid_CrossGrating,CrossGrating_Ptr,bottom,"Bottom grating","OneD_CrossGrating",0xFF);
    DEFINE_PARAMETER(Overlaid_CrossGrating,double,overlay1,"Overlay along 1st coordinate [um]","0",0xFF);
    DEFINE_PARAMETER(Overlaid_CrossGrating,double,overlay2,"Overlay along 2nd coordinate [um]","0",0xFF);
    DEFINE_PARAMETER(Overlaid_CrossGrating,double,separation,"Vertical separation between gratings [um]","0",0xFF);

    DEFINE_MODEL(Overlaid_1D_CrossGrating,CrossGrating,"One 1D grating on top of another at some angle");
    DEFINE_PTRPARAMETER(Overlaid_1D_CrossGrating,Grating_Ptr,top,"Top grating","Single_Line_Grating",0xFF);
    DEFINE_PTRPARAMETER(Overlaid_1D_CrossGrating,Grating_Ptr,bottom,"Bottom grating","Single_Line_Grating",0xFF);
    DEFINE_PARAMETER(Overlaid_1D_CrossGrating,double,angle,"Angle between gratings [deg]","90",0xFF);
    DEFINE_PARAMETER(Overlaid_1D_CrossGrating,double,separation,"Vertical separation between gratings [um]","0",0xFF);

    DEFINE_MODEL(Null_CrossGrating,CrossGrating,"A trivial crossed grating with no layers");
    DEFINE_PARAMETER(Null_CrossGrating,double,d1,"Lattice constant #1 [um]","0.5",0xFF);
    DEFINE_PARAMETER(Null_CrossGrating,double,d2,"Lattice constant #2 [um]","0.5",0xFF);
    DEFINE_PARAMETER(Null_CrossGrating,double,zeta,"Angle of lattice vectors from perpendicular [deg]","0",0xFF);

    DEFINE_MODEL(Sphere_CrossGrating,Gridded_CrossGrating,"Grating consisting of spheres embedded in a medium");
    DEFINE_PARAMETER(Sphere_CrossGrating,double,diameter,"Diameter of sphere [um]","0.05",0xFF);
    DEFINE_PARAMETER(Sphere_CrossGrating,double,above,"Distance of sphere from top of layer [um]","0.05",0xFF);
    DEFINE_PARAMETER(Sphere_CrossGrating,double,below,"Distance of sphere from bottom of layer [um]","0.05",0xFF);
    DEFINE_PARAMETER(Sphere_CrossGrating,int,nlevels,"Number of levels in structure","7",0xFF);
    DEFINE_PARAMETER(Sphere_CrossGrating,dielectric_function,sphere,"Sphere","(1,0)",0xFF);
    DEFINE_PARAMETER(Sphere_CrossGrating,dielectric_function,surrounding,"Medium around sphere","(1.5,0)",0xFF);

    DEFINE_MODEL(Pyramidal_Pit_CrossGrating,Gridded_CrossGrating,"Grating consisting of pyramidal pits");
    DEFINE_PARAMETER(Pyramidal_Pit_CrossGrating,double,side,"Length of side of base of pyramid [um]","0.05",0xFF);
    DEFINE_PARAMETER(Pyramidal_Pit_CrossGrating,double,depth,"Depth of pit [um]","0.05",0xFF);
    DEFINE_PARAMETER(Pyramidal_Pit_CrossGrating,int,nlevels,"Number of levels","10",0xFF);


}
