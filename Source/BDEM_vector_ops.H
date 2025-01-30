#ifndef BDEM_VEC_OPS_H
#define BDEM_VEC_OPS_H

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <Constants.H>

using namespace amrex;

// Define some useful stuff
// TODO: this should be reorganized
# define initV1MinusV2(V1, V2) { V1[0] - V2[0], V1[1] - V2[1], V1[2] - V2[2]  } 
# define initV1PlusV2(V1, V2) { V1[0] + V2[0], V1[1] + V2[1], V1[2] + V2[2]  } 
# define initV1PlusaByV2(V1, a, V2) { V1[0] + a*V2[0], V1[1] + a*V2[1], V1[2] + a*V2[2]  } 
# define subtractV2FromV1(V1,V2) { V1[0] -= V2[0]; V1[1] -= V2[1]; V1[2] -= V2[2]; }
# define copyVector(V1,V2) { V1[0] = V2[0]; V1[1] = V2[1]; V1[2] = V2[2]; }

AMREX_GPU_HOST_DEVICE AMREX_INLINE Real DoubleDotProd(Real T1[9],Real T2[9])
{
    Real result = 0.;

    for (int i = 0; i < 9; i++)
    {
        result += T1[i]*T2[i];
    }       

    return result;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE Real DotProd(const Real v1[3], const Real v2[3])
{
    return(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}
AMREX_GPU_HOST_DEVICE AMREX_INLINE Real normVec(const Real v1[3])
{
    return std::sqrt(DotProd(v1,v1));
}
AMREX_GPU_HOST_DEVICE AMREX_INLINE Real side_op(Real L1[6],Real L2[6])
{
    return(    L1[0]*L2[4] 
            +  L1[1]*L2[5] 
            +  L1[2]*L2[3] 
            +  L1[3]*L2[2] 
            +  L1[4]*L2[0] 
            +  L1[5]*L2[1] );
}
AMREX_GPU_HOST_DEVICE AMREX_INLINE void getvec(Real P1[3],Real P2[3],Real v[3])
{
    v[0]=P2[0]-P1[0];
    v[1]=P2[1]-P1[1];
    v[2]=P2[2]-P1[2];
}
AMREX_GPU_HOST_DEVICE AMREX_INLINE void getunitvec(Real v[3],Real vu[3])
{
    Real vmag;
    vmag=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) + TINYVAL;
    vu[0]=v[0]/vmag;
    vu[1]=v[1]/vmag;
    vu[2]=v[2]/vmag;
}
AMREX_GPU_HOST_DEVICE AMREX_INLINE void CrossProd(Real v1[3],Real v2[3],Real v[3])
{
    v[0]=v1[1]*v2[2]-v1[2]*v2[1];
    v[1]=v1[2]*v2[0]-v1[0]*v2[2];
    v[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE void generateOrthonormalBasis(Real a[3],Real b[3],Real c[3])
{
    //- Find a vector not parallel to a
    Real e[3] = {1.,0.,0.};
    Real dot_ae_sqr = DotProd(a,e);
    dot_ae_sqr *= dot_ae_sqr;

    if ( 1. - TINYVAL < dot_ae_sqr && dot_ae_sqr > 1. + TINYVAL )
    {
        e[0] = 0.;
        e[1] = 1.;
        e[2] = 2.;
    }

    //- Form vector b orthogonal to a
    CrossProd(a,e,b);
    getunitvec(b,b); 

    //- Get c orthogonal to a and b   
    CrossProd(a,b,c);
}


#endif