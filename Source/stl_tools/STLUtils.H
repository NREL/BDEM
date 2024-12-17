#ifndef STLUTILS_H
#define STLUTILS_H

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include "BDEM_vector_ops.H"

using namespace amrex;

// This gives info on the closest point to a triangle.
// In other words, it where the closest hit was. 
enum TRI_CLOSEST_TYPE
{
    POINT_A,
    POINT_B,
    POINT_C,
    EDGE_AB,
    EDGE_AC,
    EDGE_CB,
    INTERNAL    
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE void move_this_point(amrex::Real coord[3],amrex::Real newcoord[3],amrex::Real timestep,
                     int movetype,const amrex::Real u[3],const amrex::Real center[3],amrex::Real movevel)
{
  newcoord[0]=coord[0];
  newcoord[1]=coord[1];
  newcoord[2]=coord[2];
   if(movetype==1 || movetype==3)
   {
        newcoord[0] += movevel*timestep*u[0];
        newcoord[1] += movevel*timestep*u[1];
        newcoord[2] += movevel*timestep*u[2];
  
   }
   else if(movetype==2)
   {

        // Rodrigues' rotation formula
        
        amrex::Real theta = movevel*timestep;
        amrex::Real cosTheta = cos(theta);
        amrex::Real sinTheta = sin(theta);

        amrex::Real px = coord[0] - center[0];
        amrex::Real py = coord[1] - center[1]; 
        amrex::Real pz = coord[2] - center[2]; 
        
        amrex::Real rotatedX = (cosTheta + (1. - cosTheta) * u[0] * u[0]) * px +
                      ((1. - cosTheta) * u[0] * u[1] - u[2] * sinTheta) * py +
                      ((1. - cosTheta) * u[0] * u[2] + u[1] * sinTheta) * pz;

        amrex::Real rotatedY = ((1. - cosTheta) * u[1] * u[0] + u[2] * sinTheta) * px +
                        (cosTheta + (1. - cosTheta) * u[1] * u[1]) * py +
                        ((1. - cosTheta) * u[1] * u[2] - u[0] * sinTheta) * pz;

        amrex::Real rotatedZ = ((1. - cosTheta) * u[2] * u[0] - u[1] * sinTheta) * px +
                        ((1. - cosTheta) * u[2] * u[1] + u[0] * sinTheta) * py +
                        (cosTheta + (1. - cosTheta) * u[2] * u[2]) * pz;
        
        newcoord[0] = rotatedX + center[0];
        newcoord[1] = rotatedY + center[1];
        newcoord[2] = rotatedZ + center[2];

        // int t1=(movedir+1)%3;
        // int t2=(movedir+2)%3;
        // amrex::Real rad=std::pow((newcoord[t1]-center[t1]),2.0);
        // rad+=std::pow((newcoord[t2]-center[t2]),2.0);
        // rad=sqrt(rad);
   
        // amrex::Real sin_theta=(newcoord[t2]-center[t2])/rad;
        // amrex::Real cos_theta=(coord[t1]-center[t1])/rad;
        // newcoord[t1] += -movevel*rad*sin_theta*timestep;
        // newcoord[t2] +=  movevel*rad*cos_theta*timestep;
        // newcoord[movedir] += 0.0;
   }
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE void ijkFromID(int id, int gridSize[3], int ijk[3])
{
    ijk[2] = floor(id/(gridSize[0]*gridSize[1]));
    ijk[1] = floor((id - ijk[2] * (gridSize[0] * gridSize[1])) / gridSize[0]);
    ijk[0] = floor(id - ijk[2] * (gridSize[0] * gridSize[1]) - ijk[1] * gridSize[0]);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE int IdFromIJK(const int gridSize[3],const int ijk[3])
{
    return ijk[0] + ijk[1]*gridSize[0] + ijk[2]*gridSize[0]*gridSize[1];
}


AMREX_GPU_HOST_DEVICE AMREX_INLINE void ijkFromPos(Real pos[3], const Real delta[3], const Real origin[3], int ijk[3])
{
    ijk[2] = floor((pos[2] - origin[2])/delta[2]);
    ijk[1] = floor((pos[1] - origin[1])/delta[1]);
    ijk[0] = floor((pos[0] - origin[0])/delta[0]);
}

    
AMREX_GPU_HOST_DEVICE AMREX_INLINE Real Distance2(Real P1[3],Real P2[3])
{
    return( (P1[0]-P2[0])*(P1[0]-P2[0]) +
            (P1[1]-P2[1])*(P1[1]-P2[1]) +
            (P1[2]-P2[2])*(P1[2]-P2[2]) );
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE void get_plucker_coords(Real v1[3],Real v2[3],Real L[6])
{
    L[0] = v1[0]*v2[1] - v1[1]*v2[0];
    L[1] = v1[0]*v2[2] - v1[2]*v2[0];
    L[2] = v1[0]       - v2[0];
    L[3] = v1[1]*v2[2] - v1[2]*v2[1];
    L[4] = v1[2]       - v2[2];
    L[5] = v2[1]       - v1[1];
}
AMREX_GPU_HOST_DEVICE AMREX_INLINE void side_op3(Real v1[3],Real v2[3],
        Real t1[3],Real t2[3],Real t3[3],
        Real &S1, Real &S2, Real &S3)
{

    Real L[6],e1[6],e2[6],e3[6];

    get_plucker_coords(v1,v2,L);
    get_plucker_coords(t1,t2,e1);
    get_plucker_coords(t2,t3,e2);
    get_plucker_coords(t3,t1,e3);

    S1=side_op(L,e1);
    S2=side_op(L,e2);
    S3=side_op(L,e3);
}
//get normal of triangle pointing at a test-point
AMREX_GPU_HOST_DEVICE AMREX_INLINE void tri_n(Real P1[3],Real P2[3],Real P3[3],
        Real testp[3],Real n[3])
{
    Real v1[3],v2[3],magn;
    Real centr[3],c_tp_vec[3];

    getvec(P1,P2,v1);
    getvec(P1,P3,v2);
    CrossProd(v1,v2,n);

    centr[0]=0.333333*(P1[0]+P2[0]+P3[0]);
    centr[1]=0.333333*(P1[1]+P2[1]+P3[1]);
    centr[2]=0.333333*(P1[2]+P2[2]+P3[2]);

    getvec(centr,testp,c_tp_vec);
    magn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);

    if(DotProd(c_tp_vec,n) < 0.0)
    {
        magn=-magn;
    }

    n[0]=n[0]/magn;
    n[1]=n[1]/magn;
    n[2]=n[2]/magn;
}
//get normal of triangle
AMREX_GPU_HOST_DEVICE AMREX_INLINE void tri_n(Real P1[3],Real P2[3],Real P3[3],
        Real n[3])
{
    Real v1[3],v2[3],magn;

    getvec(P1,P2,v1);
    getvec(P1,P3,v2);
    CrossProd(v1,v2,n);

    magn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);

    n[0]=n[0]/magn;
    n[1]=n[1]/magn;
    n[2]=n[2]/magn;
}
AMREX_GPU_HOST_DEVICE AMREX_INLINE Real triangle_area(Real P1[3],Real P2[3],Real P3[3])
{
    Real v1[3],v2[3],area[3];

    getvec(P1,P2,v1);
    getvec(P1,P3,v2);
    CrossProd(v1,v2,area);
    return(0.5 * sqrt(DotProd(area,area)) );
}
AMREX_GPU_HOST_DEVICE AMREX_INLINE Real triangle_area_2d(Real P1[2],Real P2[2],Real P3[2])
{
    //(x1-x2)i+(y1-y2)j x (x1-x3)i+(y1-y3)j
    //(x1-x2)*(y1-y3)+(x3-x1)*(y1-y2)
    //y3*(x2-x1)+y1*(x1-x2)+y2*(x1-x3)+y1*(x3-x1)
    //y3*(x2-x1)+y1*(x3-x2)+y2*(x1-x3)
    //-1*(y1*(x2-x3)+y2*(x3-x1)+y3*(x1-x2))
    //
    //
    Real area=0.5*(P1[1]*(P2[0]-P3[0])+P2[1]*(P3[0]-P1[0])+P3[1]*(P1[0]-P2[0]));
    return(amrex::Math::abs(area));
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE Real distsq_2d(Real P1[2],Real P2[2])
{
    return( (P1[0]-P2[0])*(P1[0]-P2[0]) + (P1[1]-P2[1])*(P1[1]-P2[1]) );
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE Real distsq_from_line_2d(Real P0[2],Real P1[2],Real P2[2])
{
    Real dist,dist2;

    dist = (P2[1]-P1[1])*(P0[0]-P1[0]) - (P2[0]-P1[0])*(P0[1]-P1[1]);
    dist = dist/sqrt( (P2[1]-P1[1])*(P2[1]-P1[1]) + (P2[0]-P1[0])*(P2[0]-P1[0]) );
    dist2=dist*dist;

    return(dist2);
}
//this is only useful when v1-v2 segment intersects the triangle
AMREX_GPU_HOST_DEVICE AMREX_INLINE bool find_intersection_point(Real v1[3],Real v2[3],
        Real t1[3], Real t2[3], Real t3[3],Real ip[3],int bisect_iters=20,Real tol=1e-6)
{
    Real plane_eq_mid,plane_eq1,plane_eq2;

    Real ab[3],ac[3],n[3],magn;
    Real midp[3],p1[3],p2[3];

    getvec(t1,t2,ab);
    getvec(t1,t3,ac);

    CrossProd(ab,ac,n);
    magn=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);

    n[0]=n[0]/magn;
    n[1]=n[1]/magn;
    n[2]=n[2]/magn;

    p1[0]=v1[0];
    p1[1]=v1[1];
    p1[2]=v1[2];

    p2[0]=v2[0];
    p2[1]=v2[1];
    p2[2]=v2[2];


    bool all_ok=true;

    for(int i=0;i<bisect_iters;i++)
    {
        midp[0]=0.5*(p1[0]+p2[0]);
        midp[1]=0.5*(p1[1]+p2[1]);
        midp[2]=0.5*(p1[2]+p2[2]);

        plane_eq_mid= (midp[0]-t1[0])*n[0] + (midp[1]-t1[1])*n[1] + (midp[2]-t1[2])*n[2];
        plane_eq1   = (p1[0]  -t1[0])*n[0] + (p1[1]  -t1[1])*n[1] + (p1[2]  -t1[2])*n[2];
        plane_eq2   = (p2[0]  -t1[0])*n[0] + (p2[1]  -t1[1])*n[1] + (p2[2]  -t1[2])*n[2];

        //Print()<<"midp:"<<midp[0]<<"\t"<<midp[1]<<"\t"<<midp[2]<<"\t"<<plane_eq_mid<<"\n";

        if(fabs(plane_eq_mid) < tol)
        {
            break;
        }

        if(plane_eq_mid*plane_eq1 < 0.0)
        {
            p2[0]=midp[0];
            p2[1]=midp[1];
            p2[2]=midp[2];
        }
        else if(plane_eq_mid*plane_eq2 < 0.0)
        {
            p1[0]=midp[0];
            p1[1]=midp[1];
            p1[2]=midp[2];
        }
        else //plane_eq_mid is 0 
            //or error: p1,midp and p2 are on the same side 
            //which is not what this function is meant for
        {
            if(plane_eq_mid*plane_eq1 > 0.0 and plane_eq_mid*plane_eq2 > 0.0)
            {
                all_ok=false;
            }
            break;
        }
    }

    ip[0]=midp[0];
    ip[1]=midp[1];
    ip[2]=midp[2];

    return(all_ok);
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE int lineseg_tri_intersect(Real v1[3],Real v2[3],
        Real t1[3],Real t2[3],Real t3[3])
{
    //see plucker coordinates based method
    //https://members.loria.fr/SLazard/ARC-Visi3D/Pant-project/files/Line_Triangle.html

    Real S1,S2,S3;
    Real tri_area,area1,area2;
    Real L2[6],L3[6],L4[6],ls_s1,ls_s2;

    side_op3(v1,v2,t1,t2,t3,S1,S2,S3);

    //we are assuming there are no intersections initially
    int no_intersections=1;

    Real eps = std::numeric_limits<Real>::epsilon();

    //coplanar (S1,S2,S3 = 0)
    if(fabs(S1) < eps and fabs(S2) < eps and fabs(S3) < eps) 
    {
        //Print()<<"line segment and triangle are in the same plane\t"<<S1<<"\t"<<S2<<"\t"<<S3<<"\n";
        tri_area=triangle_area(t1,t2,t3);

        /*if(tri_area == 0)
            {
            amrex::Abort("problem with triangle\n");
            }*/
        area1=(triangle_area(t1,t2,v1)+triangle_area(t2,t3,v1)+triangle_area(t3,t1,v1));
        area2=(triangle_area(t1,t2,v2)+triangle_area(t2,t3,v2)+triangle_area(t3,t1,v2));

        if( fabs(area1-tri_area)>eps or fabs(area2-tri_area)>eps )
        {
            no_intersections = 0;
        }
    }
    //proper and edge intersection
    else if( (S1 < 0.0 and S2 < 0.0 and S3 < 0.0) or 
            (S1 > 0.0 and S2 > 0.0 and S3 > 0.0) or
            (fabs(S1) < eps and S2*S3 > 0.0) or     //S1=0
            (fabs(S2) < eps and S3*S1 > 0.0) or     //S2=0
            (fabs(S3) < eps and S1*S2 > 0.0) )      //S3=0
    {

        get_plucker_coords(v1,t1,L2);
        get_plucker_coords(t1,v2,L3);
        get_plucker_coords(t2,t3,L4);

        /*if(fabs(S1*S2*S3) < eps)
            {
            Print()<<"edge intersection S1,S2,S3:"
            <<S1<<"\t"<<S2<<"\t"<<S3<<"\n";
            }*/

        ls_s1 = side_op(L4,L3);
        ls_s2 = side_op(L4,L2);

        if(ls_s1*ls_s2 > 0.0)
        {
            no_intersections = 0;
        }
    }
    //vertex intersection
    else if( (fabs(S1) < eps and fabs(S2) < eps) or  //S1,S2=0
            (fabs(S2) < eps and fabs(S3) < eps) )   //S2,S3=0
    {

        //Print()<<"vertex intersection type 1\n";   
        //don't chose vertex 2 or 3
        get_plucker_coords(v2,t1,L2);
        get_plucker_coords(t1,v1,L3);
        get_plucker_coords(t2,t3,L4);

        ls_s1=side_op(L4,L3);
        ls_s2=side_op(L4,L2);

        if(ls_s1*ls_s2 > 0)
        {
            no_intersections = 0;
        }
    }
    else if(fabs(S3) < eps and fabs(S1) < eps) //S3,S1=0
    {

        //Print()<<"vertex intersection type 2\n";
        //don't chose vertex 1
        get_plucker_coords(v2,t2,L2);
        get_plucker_coords(t2,v1,L3);
        get_plucker_coords(t3,t1,L4);

        ls_s1=side_op(L4,L3);
        ls_s2=side_op(L4,L2);

        if(ls_s1*ls_s2 > 0)
        {
            no_intersections=0;
        }
    }

    return(no_intersections);

}
AMREX_GPU_HOST_DEVICE AMREX_INLINE Real point_tri_distance(Real p0[3],
        Real rot_pts[9],int domndir,Real rotmat[9])
{
    Real p0t[3],p1t[3],p2t[3],p3t[3];
    Real t0[2];
    Real t1[2];
    Real t2[2];
    Real t3[2];
    Real tol=1e-9;
    Real ptdist=0.0;
    Real mindist1sq,mindist2sq,mindist3sq,mindistsq,mindistsq_l;
    Real v01sq,v02sq,v03sq;
    Real s1,s2,s3;
    Real s01sq,s02sq,s03sq;
    int sign_inside,sign_outside;
    int sign_t1t0,sign_t3t0,sign_t1t3,sign_t3t1;
    int sign_t2t0,sign_t2t3,sign_t3t2;
    Real bigval=1e50;

    int dir1,dir2;

    dir1=(domndir+1)%3;
    dir2=(domndir+2)%3;

    //rotation
    for(int i=0;i<3;i++)
    {
        p0t[i]=0.0;
        for(int j=0;j<3;j++)
        {
            p0t[i] += rotmat[i*3+j]*p0[j];
        }
    }
    for(int i=0;i<3;i++)
    {
        p1t[i]=rot_pts[i];
        p2t[i]=rot_pts[3+i];
        p3t[i]=rot_pts[6+i];
    }

    //translation
    for(int i=0;i<3;i++)
    {
        p0t[i] -= p1t[i];
        p2t[i] -= p1t[i];
        p3t[i] -= p1t[i];

        p1t[i] -= p1t[i];
    }

    t0[0]=p0t[dir1];
    t0[1]=p0t[dir2];

    t1[0]=p1t[dir1];
    t1[1]=p1t[dir2];

    t2[0]=p2t[dir1];
    t2[1]=p2t[dir2];

    t3[0]=p3t[dir1];
    t3[1]=p3t[dir2];

    //find barycentric coords
    //now t1 is (0,0), t2 (a,0) and t3 (b,c)
    //
    Real a=t2[0];
    Real b=t3[0];
    Real c=t3[1];
    //t0=r1*t1+r2*t2+r3*t3
    //(r1+r2+r3)=1

    Real r3=t0[1]/c;
    Real r2=(t0[0]-r3*b)/a;
    Real r1=1.0-r2-r3;

    if(r1*(1-r1)>0.0 && r2*(1-r2)>0.0 && r3*(1-r3)>0.0)
    {
        ptdist=amrex::Math::abs(p0t[domndir]);
    }
    else
    {
        //check t1t2 line
        sign_inside=amrex::Math::copysign(1.0,t3[1]);
        sign_outside=amrex::Math::copysign(1.0,t0[1]);

        if(sign_inside*sign_outside<0.0 && t0[0]*(t0[0]-a) < tol)
        {
            mindistsq_l=t0[1]*t0[1];
        }
        else
        {
            //Print()<<"checking t1t3\n";
            sign_inside=amrex::Math::copysign(1.0,-c*a);
            sign_outside=amrex::Math::copysign(1.0,b*t0[1]-c*t0[0]);

            sign_t1t0=amrex::Math::copysign(1.0,c*t0[1] + b*t0[0]);
            sign_t1t3=amrex::Math::copysign(1.0,c*t3[1] + b*t3[0]);
            sign_t3t0=amrex::Math::copysign(1.0,c*(t0[1]-c) + b*(t0[0]-b));
            sign_t3t1=amrex::Math::copysign(1.0,c*(t1[1]-c) + b*(t1[0]-b));

            if(sign_inside*sign_outside<=0.0 && sign_t1t0*sign_t1t3>=0.0 && sign_t3t0*sign_t3t1>=0.0)
            {
                mindistsq_l=distsq_from_line_2d(t0,t3,t1);
            }
            else
            {
                //Print()<<"checking t2t3\n";
                sign_inside=amrex::Math::copysign(1.0,c*a);
                sign_outside=amrex::Math::copysign(1.0,(b-a)*t0[1]-c*(t0[0]-a));

                sign_t2t0=amrex::Math::copysign(1.0,c*t0[1] + (b-a)*(t0[0]-a));
                sign_t2t3=amrex::Math::copysign(1.0,c*t3[1] + (b-a)*(t3[0]-a));
                sign_t3t0=amrex::Math::copysign(1.0,c*(t0[1]-c) + (b-a)*(t0[0]-b));
                sign_t3t2=amrex::Math::copysign(1.0,c*(t2[1]-c) + (b-a)*(t2[0]-b));

                if(sign_inside*sign_outside<=0.0 && sign_t2t0*sign_t2t3>=0.0 && sign_t3t0*sign_t3t2>=0.0)
                {
                    mindistsq_l=distsq_from_line_2d(t0,t2,t3);
                }
                else
                {
                    mindistsq_l=bigval;
                }
            }
        }

        //vertex distances
        v01sq = t0[0]*t0[0] + t0[1]*t0[1];
        v02sq = (t0[0]-a)*(t0[0]-a) + t0[1]*t0[1];
        v03sq = (t0[0]-b)*(t0[0]-b) + (t0[1]-c)*(t0[1]-c);

        mindistsq=amrex::min(mindistsq_l,amrex::min(v01sq,amrex::min(v02sq,v03sq)));

        ptdist=sqrt(p0t[domndir]*p0t[domndir] + mindistsq);
    }

    return(ptdist);
} 


AMREX_GPU_HOST_DEVICE AMREX_INLINE Real point_tri_distance(Real P[3],
        Real A[3],Real B[3], Real C[3], int& closestType, Real dist[3] )
{

    // Adapted from
    // Real-time collision detection, Christer Ericson, 2005, p136-142

    // This is an algorithm developed for videogames and adapted from
    // OpenFOAM-v2112

    //* How does it work?
    //* The method subdivides the triangle into seven sub-regions
    //*       \| B
    //*        |\
    //*        | \
    //*        |  \
    //*      __|_ _\__
    //*       A|   C\
    //* Three sub-regions are outside vertices (vertex is the closest point)
    //* Three sub-regions are outside edges (point on edge is closest)
    //* One sub-region is inside the triangle  (point given by barycentric is closest)

    // Check if the projection of P is in the region outside A
        Real AB[3] = initV1MinusV2(B,A);
        Real AC[3] = initV1MinusV2(C,A);
        Real AP[3] = initV1MinusV2(P,A);

        // These are projections over the edges
        Real d1 = DotProd(AB,AP);
        Real d2 = DotProd(AC,AP);

        if ( d1 <= 0. && d2 <= 0. )
        {
            // The closest point is A
            // Barycentric coordinates (1 0 0)
            copyVector(dist,AP)
            closestType = TRI_CLOSEST_TYPE::POINT_A;
            return normVec(AP);
        }

    // Check if the projection of P is in region outside B
        Real BP[3] = initV1MinusV2(P,B);
        Real d3 = DotProd(AB,BP);
        Real d4 = DotProd(AC,BP);

        if ( d3 >= 0. && d4 <= d3 )
        {
            // The closest point is B
            // Barycentric coordinates (0 1 0)
            copyVector(dist,BP)
            closestType = TRI_CLOSEST_TYPE::POINT_B;
            return normVec(BP);
        }
    
    // Check if the projection of P is in the region outside AB
        Real vc = d1*d4 - d3*d2;

        if ( vc <= 0. && d1 >= 0. && d3 <= 0. )
        {
            //! Assume non-degenerate triangle
            //! You get printstack if the triangle is degenerate

            // Barycentric coordinates (1-v,v,0)
            Real v = d1 / ( d1 - d3 );

            Real nearP[3] = initV1PlusaByV2(A,v,AB);
            Real nearDist[3] = initV1MinusV2(P,nearP);

            copyVector(dist,nearDist)
            closestType = TRI_CLOSEST_TYPE::EDGE_AB;
            return normVec(dist);
        }

    // Check if the projection of P is in region outside C
        Real CP[3] = initV1MinusV2(P,C);
        Real d5 = DotProd(AB,CP);
        Real d6 = DotProd(AC,CP);

        if ( d6 >= 0. && d5 <= d6 )
        {
            // The closest point is C
            // Barycentric coordinates (0 1 0)
            copyVector(dist,CP)
            closestType = TRI_CLOSEST_TYPE::POINT_C;
            return normVec(CP);
        }

    // Check if the projection of P is in the region outside AC
        Real vb = d5*d2 - d1*d6;

        if ( vb <= 0. && d2 >= 0. && d6 <= 0. )
        {
            //! Assume non-degenerate triangle
            //! You get printstack if the triangle is degenerate

            // Barycentric coordinates (1-w,0,w)
            Real w = d2 / ( d2 - d6 );

            Real nearP[3] = initV1PlusaByV2(A,w,AC);
            Real nearDist[3] = initV1MinusV2(P,nearP);
            copyVector(dist,nearDist)
            closestType = TRI_CLOSEST_TYPE::EDGE_AC;
            return normVec(dist);
        }

    // Check if the projection of P is in the region outside BC
        Real va = d3*d6 - d5*d4;

        if ( va <= 0. && d4 >= d3 && d5 >= d6 )
        {
            //! Assume non-degenerate triangle
            //! You get printstack if the triangle is degenerate

            // Barycentric coordinates (0,1-w,w)
            Real w = ( d4 - d3 ) / ( ( d4 - d3 ) + ( d5 - d6 ));
            Real CB[3] = initV1MinusV2(C,B);
            Real nearP[3] = initV1PlusaByV2(B,w,CB);
            Real nearDist[3] = initV1MinusV2(P,nearP);
            copyVector(dist,nearDist)
            closestType = TRI_CLOSEST_TYPE::EDGE_CB;
            return normVec(dist);
        }
    
    // If inside the triangle, use barycentric coordinates
        Real denom = 1. / ( va + vb + vc );
        Real v = vb * denom;
        Real w = vc * denom;

        Real nearDist[3];
        for (int i = 0; i < 3; i++)
        {
            nearDist[i] = P[i] - ( A[i] + v * AB[i] + w * AC[i] );  
        }
        copyVector(dist,nearDist)
        closestType = TRI_CLOSEST_TYPE::INTERNAL;
        return normVec(dist);
        
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE Real point_tri_distance(Real P[3],
        Real AA[3],Real BB[3], Real CC[3], Real PCP[3])
{
    Real E0[3],E1[3],DD[3], CP[3];
    Real a,b,c;
    Real dist2;
    Real n[3];
    Real s;
    Real t;

    a=b=c=0.0;
    for (int d=0;d<3;d++) 
    {
        E0[d]=CC[d]-BB[d];
        E1[d]=AA[d]-BB[d];
        a+=(E0[d]*E0[d]);
        b+=(E0[d]*E1[d]);
        c+=(E1[d]*E1[d]);   
        DD[d]=BB[d]-P[d];
    }  
    n[0]=E0[1]*E1[2]-E0[2]*E1[1];
    n[1]=E0[2]*E1[0]-E0[0]*E1[2];
    n[2]=E0[0]*E1[1]-E0[1]*E1[0];
    Real det=a*c-b*b;
    Real detinv=1.0/det;
    Real denom=a-2*b+c;
    Real dp=DD[0]*n[0]+DD[1]*n[1]+DD[2]*n[2];
    Real d=E0[0]*DD[0]+E0[1]*DD[1]+E0[2]*DD[2];
    Real e=E1[0]*DD[0]+E1[1]*DD[1]+E1[2]*DD[2];
    s=b*e-c*d;
    t=b*d-a*e;
    int iregion;  
    if ((s+t) <= det) {
        if (s < 0) {
            if (t < 0) {
                iregion=4;
            } else {
                iregion=3;
            }
        }
        else if (t < 0) {
            iregion=5;
        } else {
            iregion=0;
        }
    } else {
        if (s < 0) {
            iregion=2;
        } else if (t < 0) {
            iregion=6;
        } else {
            iregion=1;
        }
    }

    Real numer,tmp0,tmp1;

    switch(iregion) {
        case 0:
            s=s*detinv;
            t=t*detinv;
            break;
        case 1:
            numer=c+e-b-d;
            if (numer <=0) {
                s=0;
            } else {
                if (numer >=denom) {
                    s=1;
                } else {
                    s=numer/denom;
                }
            }
            t=1-s;
            break;
        case 2:
            tmp0=b+d;
            tmp1=c+e;
            if (tmp1 > tmp0) {
                numer=tmp1-tmp0;
                if (numer >=denom) {
                    s=1;
                } else {
                    s=numer/denom;
                }
                t=1-s;
            }
            else  {
                s=0;
                if (tmp1 <=0) {
                    t=1;
                } else if (e >=0) {
                    t=0;
                } else {
                    t=-e/c;
                }
            }
            break;
        case 3:
            s=0;
            if (e>=0) {
                t=0;
            } else if (c+e <=0) {
                t=1;
            } else {
                t=-e/c;
            }
            break;
        case 4:
            if (d < 0) {
                t=0;
                if (a+d <=0) {
                    s=1;
                } else {
                    s=-d/a;
                }
            } else {
                s=0;
                if (c+e <=0) {
                    t=1;
                } else if (e>=0) {
                    t=0;
                } else {
                    t=-e/c;
                }
            }
            break;
        case 5:
            t=0;
            if (d>=0) {
                s=0;
            } else if (a+d <=0) {
                s=1;
            } else {
                s=-d/a;
            }
            break;
        case 6:
            tmp0=b+e;
            tmp1=a+d;
            if (tmp1>tmp0) {
                numer=tmp1-tmp0;
                if (numer>=denom) {
                    t=1;
                } else {
                    t=numer/denom;
                }
                s=1-t;
            } else {
                t=0;
                if (tmp1<=0) {
                    s=1;
                } else if (d>=0) {
                    s=0;
                } else {
                    s=-d/a;
                }
            }
    }
    dist2=0.0;
    for(int d=0;d<3;d++) 
    {
        CP[d]=BB[d]+s*E0[d]+t*E1[d]; 
        PCP[d] = P[d] - CP[d];   
        dist2+=(CP[d]-P[d])*(CP[d]-P[d]);
    }

    Real dist=sqrt(dist2);

    return(dist);
}



AMREX_GPU_HOST_DEVICE AMREX_INLINE Real point_tri_distance(Real P[3],
        Real AA[3],Real BB[3], Real CC[3])
{
    Real CP[3];
    return point_tri_distance(P, AA, BB, CC, CP);
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void qrfac(double A[9],double Q[9])
{
    for(int i=0;i<3;i++)
    {
        for(int dim=0;dim<3;dim++)
        {
            Q[3*i+dim]=A[3*i+dim];
        }
        for(int v=0;v<i;v++)
        {
            double dotpdt=0.0;
            for(int dim=0;dim<3;dim++)
            {
                dotpdt += A[3*i+dim]*Q[3*v+dim];
            }
            for(int dim=0;dim<3;dim++)
            {
                Q[3*i+dim] -= dotpdt*Q[3*v+dim];
            }
        }

        //find magnitude
        double mag=0.0;
        for(int dim=0;dim<3;dim++)
        {
            mag += Q[3*i+dim]*Q[3*i+dim];
        }
        mag=sqrt(mag);


        for(int dim=0;dim<3;dim++)
        {
            Q[3*i+dim]/=mag;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void matmult(double A[9],double B[9],double C[9])
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            C[3*j+i]=0.0;
            for(int k=0;k<3;k++)
            {
                C[3*j+i] += A[3*k+i]*B[3*j+k];
            }
        }
    }
}
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void transpose(double A[9],double At[9])
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            At[3*i+j]=A[3*j+i];
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void get_qtaq(double A[9],double Q[9],double qtaq[9])
{
    double temp[9],Qt[9];
    matmult(A,Q,temp);
    transpose(Q,Qt);
    matmult(Qt,temp,qtaq);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void geteigvectors(double A[9],double eigmat[9])
{
    double Anext[9]={0.0};
    double Q[9]={0.0};
    double prod[9]={1,0,0,0,1,0,0,0,1};
    int maxiter=100;

    for(int i=0;i<9;i++)
        eigmat[i]=prod[i];

    for(int iter=0;iter<maxiter;iter++)
    {
        qrfac(A,Q);
        get_qtaq(A,Q,Anext);
        matmult(eigmat,Q,prod);
        for(int i=0;i<9;i++)
        {
            A[i]=Anext[i];
            eigmat[i]=prod[i];
        }
    }
}


#endif