//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: brdf.cpp
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
#include <complex>
#include <limits>
#include "scatmech.h"
#include "brdf.h"
#include "askuser.h"

using namespace std;


namespace SCATMECH {


    //
    // The following function gives the rotation angles for the
    // incident and scattered polarization coordinates to change from
    // a {s,p} coordinate system on the incident and scattered light
    // to a coordinate system whereby the polarization coordinates are
    // parallel and perpendicular to the plane subtended by the incident
    // and scattered directions.
    //
    // This function was added in Version 3 (TAG: 18 JAN 2001)
    //
    void
    Coordinate_System_plane_angles(double &psii,double &psis,
                                   double thetai,double thetas,double phis)
    {
        double costhetai=cos(thetai);
        double costhetas=cos(thetas);
        double cosphi=cos(phis);
        double sinthetai=sin(thetai);
        double sinthetas=sin(thetas);
        double sinphi=sin(phis);

        double cosd = costhetai*costhetas - cosphi*sinthetai*sinthetas;
        double sind = sqrt(1.-sqr(cosd));
        double tana = (costhetas-costhetai*cosd)/(costhetai*sind);
        double sina = tana/sqrt(1. + sqr(tana));
        double cosa = 1./sqrt(1. + sqr(tana));
        double sinb = -sinphi*sinthetas*sinthetai/sind;
        double cosb = sqrt(1. - sqr(sinb));
        double sinad = cosd*sina-cosa*sind;

        psii = -atan2(-sinb,cosb*sina);
        psis = atan2(sinb,cosb*sinad);
    }

    int
    BRDF_Model::
    set_geometry(double _thetai,double _thetas,double _phis,double _rotation)
    {
        thetai=_thetai;
        thetas=_thetas;
        phis=_phis;
        rotation=_rotation;

        if (thetai!=last_thetai || rotation!=last_rotation) set_recalc(recalc_on_incident_change);
        if (thetas!=last_thetas || phis!=last_phis || rotation!=last_rotation) set_recalc(recalc_on_scatter_change);

        last_thetai=thetai;
        last_thetas=thetas;
        last_phis=phis;
        last_rotation=rotation;

        if (thetas>=pi/2) return 0;
        if (thetas<=-pi/2) return 0;
        if (thetai>=pi/2) return 0;
        if (thetai<=-pi/2) return 0;
        return 1;
    }

    int
    BRDF_Model::
    set_geometry(const Vector& source,const Vector& viewer,const Vector& normal,const Vector& xaxis)
    {
        // The vectors need not be normalized...
        Vector norm_source = unit(source);
        Vector norm_viewer = unit(viewer);
        Vector norm_normal = unit(normal);
        // The x-axis does not need to be orthogonal to normal
        Vector norm_xaxis = unit(xaxis-norm_normal*(norm_normal*xaxis));

        // cosine of the incident angle
        double costhetai = norm_source*norm_normal;
        // cosine of the scattering angle
        double costhetas = norm_viewer*norm_normal;

        // Get scattering and incident angles...
        if (((costhetai<0. || costhetas<0.) && type==Type_DOWNUP) ||
                ((costhetai<0. || costhetas>0.) && type==Type_DOWNDOWN) ||
                ((costhetai>0. || costhetas>0.) && type==Type_UPDOWN) ||
                ((costhetai>0. || costhetas<0.) && type==Type_UPUP)) {

            thetas = thetai = 0.;
            phis = 0.;

            if (thetai!=last_thetai || rotation!=last_rotation) set_recalc(recalc_on_incident_change);
            if (thetas!=last_thetas || phis!=last_phis || rotation!=last_rotation) set_recalc(recalc_on_scatter_change);

            last_thetai=thetai;
            last_thetas=thetas;
            last_phis=phis;
            last_rotation=rotation;

            return (0);
        }

        if (costhetai>1.0) costhetai=1.;
        if (costhetas>1.0) costhetas=1.;

        thetai = acos(costhetai);
        thetas = acos(costhetas);

        // Project kin and kout onto surface plane
        Vector kinperp  = norm_source - costhetai*norm_normal;
        Vector koutperp = norm_viewer - costhetas*norm_normal;

        // Make them unit vectors...
        Vector kinperphat = unit(kinperp);

        // Get a vector perpendicular to kinperphat and normal...
        Vector yhat = perpto(kinperphat,norm_normal);

        double cosphis = -norm_viewer*kinperphat;
        double sinphis = -norm_viewer*yhat;

        phis = atan2(sinphis,cosphis);

        Vector norm_yaxis = perpto(norm_normal,norm_xaxis);

        rotation = atan2(-norm_yaxis*norm_source,norm_xaxis*norm_source);

        if (thetai!=last_thetai || rotation!=last_rotation) set_recalc(recalc_on_incident_change);
        if (thetas!=last_thetas || phis!=last_phis || rotation!=last_rotation) set_recalc(recalc_on_scatter_change);

        last_thetai=thetai;
        last_thetas=thetas;
        last_phis=phis;
        last_rotation=rotation;

        return 1;
    }

    void
    BRDF_Model::
    convert(JonesMatrix& J,Coordinate_System to) const
    {
        Coordinate_System from = model_cs;
        if (from==to) return;

        if (from==xyxy) {
            // Convert to BRDF_Model::psps
            J = JonesRotator(phis)*J;
        } else if (from==plane) {
            // Convert to BRDF_Model::psps
            double psii,psis;
            Coordinate_System_plane_angles(psii,psis,thetai,thetas,phis);
            J = (JonesRotator(psis)*J)*JonesRotator(psii);
        }

        if (to == xyxy) {
            // Convert to BRDF_Model::xyxy
            J = JonesRotator(-phis)*J;
        } else if (to==plane) {
            // Convert to BRDF_Model::plane
            double psii,psis;
            Coordinate_System_plane_angles(psii,psis,thetai,thetas,phis);
            J = (JonesRotator(-psis)*J)*JonesRotator(-psii);
        }
    }

    void
    BRDF_Model::
    convert(MuellerMatrix& M,Coordinate_System to) const
    {
        Coordinate_System from = model_cs;
        if (from==to) return;

        if (from==xyxy) {
            // Convert to BRDF_Model_psps
            M = ((MuellerMatrix)JonesRotator(phis))*M;
        } else if (from==plane) {
            // Convert to BRDF_Model_psps
            double psii,psis;
            Coordinate_System_plane_angles(psii,psis,thetai,thetas,phis);
            M= (((MuellerMatrix)JonesRotator(psis))*M)*
               ((MuellerMatrix)JonesRotator(psii));
        }

        if (to == xyxy) {
            // Convert to BRDF_Model_xyxy
            M = ((MuellerMatrix)JonesRotator(-phis))*M;
        } else if (to==plane) {
            // Convert to BRDF_Model_plane
            double psii,psis;
            Coordinate_System_plane_angles(psii,psis,thetai,thetas,phis);
            M= (((MuellerMatrix)JonesRotator(-psis))*M)*
               ((MuellerMatrix)JonesRotator(-psii));
        }
    }


    //
    // The following function returns the jones matrix for scattering using a
    // specific coordinate system, presumably, but not necessarily different,
    // than the {s,p} system.
    //
    // This function was added in Version 3 (TAG: 18 JAN 2001)
    //
    JonesMatrix
    BRDF_Model::
    Jones(double thetai,double thetas,double phis,double rotation,Coordinate_System cs)
    {
        if (set_geometry(thetai,thetas,phis,rotation)) {
            JonesMatrix j = jones();
            convert(j,cs);
            return j;
        } else {
            return JonesZero();
        }
    }

    //
    // The following function returns the mueller matrix for scattering using a
    // specific coordinate system, presumably, but not necessarily different,
    // than the {s,p} system.
    //
    // This function was added in Version 3 (TAG: 18 JAN 2001)
    //
    MuellerMatrix
    BRDF_Model::
    Mueller(double thetai,double thetas,double phis,double rotation,Coordinate_System cs)
    {
        if (set_geometry(thetai,thetas,phis,rotation)) {
            MuellerMatrix m = mueller();
            convert(m,cs);
            return m;
        } else {
            return MuellerZero();
        }
    }

    void
    BRDF_Model::
    setup()
    {
        Model::setup();
        if (type!=Type_DOWNUP && substrate.k(lambda)!=0.)
            error("Absorbing substrate for type==1, 2, or 3");
        if (type<0||type>3)
            error("Illegal type: must be 0, 1, 2, or 3");
    }

    JonesMatrix
    BRDF_Model::
    Jones(const Vector& source,
          const Vector& viewer,
          const Vector& normal,
          const Vector& xaxis,
          Coordinate_System cs)
    {
        if (set_geometry(source,viewer,normal,xaxis)) {
            return Jones(thetai,thetas,phis,rotation,cs);
        } else {
            return JonesZero();
        }
    }

    MuellerMatrix
    BRDF_Model::
    Mueller(const Vector& source,
            const Vector& viewer,
            const Vector& normal,
            const Vector& xaxis,
            Coordinate_System cs)
    {
        if (set_geometry(source,viewer,normal,xaxis)) {
            return Mueller(thetai,thetas,phis,rotation,cs);
        } else {
            return MuellerZero();
        }
    }

    DEFINE_VIRTUAL_MODEL(BRDF_Model,Model,"All BRDF models.");
    DEFINE_PARAMETER(BRDF_Model,double,lambda,"Wavelength [um]","0.532",0xFF);
    DEFINE_PARAMETER(BRDF_Model,dielectric_function,substrate,"Substrate","(4.05,0.05)",0xFF);
    DEFINE_PARAMETER(BRDF_Model,int,type,"(0) for Forward/Reflection, (1) for Forward/Transmission, (2) for Backward/Reflection, or (3) for Backward/Transmission","0",0xFF);

    int get_model_depth=0;


} // namespace SCATMECH


