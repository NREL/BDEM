/* ============================================================================
    This file is part of BDEM.

    Description: 
    Class for managing one body described by a STL file.
============================================================================ */

#ifndef __BDEM_STL_BODY_H__
#define __BDEM_STL_BODY_H__

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>

namespace BDEM
{

class STLBody
{
    private:

        /* Private members */
            const amrex::ParmParse&                     parameters;             //- Main input file
            std::string                                 stl_name;               //- Name of the STL

            AMREX_GPU_DEVICE_MANAGED int                num_tri;                //- Number of trinagles 
            AMREX_GPU_DEVICE_MANAGED int                ndata_per_tri;
            AMREX_GPU_DEVICE_MANAGED int                ndata_per_normal;
            AMREX_GPU_DEVICE_MANAGED amrex::Real        eigdirs[9];
            AMREX_GPU_DEVICE_MANAGED amrex::Real        bbox_lo[3];
            AMREX_GPU_DEVICE_MANAGED amrex::Real        bbox_hi[3];
            AMREX_GPU_DEVICE_MANAGED amrex::Real        expfac;
            int                                         nlines_per_facet;

            /* Array of triangle point coordinates (x,y,z) */
            AMREX_GPU_DEVICE_MANAGED amrex::Real*       tri_pts; 
            
            /* Array of triangle coord (AoS) */
            amrex::Gpu::ManagedVector<amrex::Real>*     tri_pts_vec;

            /* Array of triangle normals (x,y,z)*/
            AMREX_GPU_DEVICE_MANAGED amrex::Real*       tri_normals; 
            
            /* Array of triangle point coordinates (AoS)*/
            amrex::Gpu::ManagedVector<amrex::Real>*     tri_normals_vec;


            // TODO: I need to understand what the rest of this is
            AMREX_GPU_DEVICE_MANAGED amrex::Real*       tri_mat;
            
            amrex::Gpu::ManagedVector<amrex::Real>*     tri_mat_vec; 
            
            AMREX_GPU_DEVICE_MANAGED int*               domndir;
            
            amrex::Gpu::ManagedVector<int>*             domndir_vec;

            AMREX_GPU_DEVICE_MANAGED amrex::Real*       rot_tri_pts;
            
            amrex::Gpu::ManagedVector<amrex::Real>*     rot_tri_pts_vec; 

            AMREX_GPU_DEVICE_MANAGED int*               sorted_indexarray;

            amrex::Gpu::ManagedVector<int>*             sorted_indexarray_vec; 

            std::map < 
                std::pair<int,int>, 
                amrex::Array<amrex::Real,6> 
                >                                       boxmap;
        

        /* Private member functions */

            /* Read STL from file */
            void readFromFile();

    public:

        /* Constructors and destructor */
        
            STLBody() = delete;     // Disallow default constructor

            /* Construct by input file and name of the STL object. */
            STLBody( const amrex::ParmParse& input, std::string stlname );

            ~STLBody();

        /* Public member functions */

            /* Write to new STL file. Takes the time step as argument. */
            void write( int timeStep ) const;

};

}

#endif // MACRO
