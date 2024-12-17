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

using namespace amrex;

// Define some useful stuff
// TODO: this should be reorganized
# define initV1MinusV2(V1, V2) { V1[0] - V2[0], V1[1] - V2[1], V1[2] - V2[2]  } 
# define initV1PlusaByV2(V1, a, V2) { V1[0] + a*V2[0], V1[1] + a*V2[1], V1[2] + a*V2[2]  } 
# define subtractV2FromV1(V1,V2) { V1[0] -= V2[0]; V1[1] -= V2[1]; V1[2] -= V2[2]; }
# define copyVector(V1,V2) { V1[0] = V2[0]; V1[1] = V2[1]; V1[2] = V2[2]; }

AMREX_GPU_HOST_DEVICE AMREX_INLINE Real DotProd(Real v1[3],Real v2[3])
{
    return(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}
AMREX_GPU_HOST_DEVICE AMREX_INLINE Real normVec(Real v1[3])
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
    vmag=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
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

#endif