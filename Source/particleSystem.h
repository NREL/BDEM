/**
 * Author: Federico Municchi
 * 
 * Description: 
 * Base class for managing particles including particle movement.
 */

#pragma once

#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>
#include <AMReX_ParmParse.H>


//- Keep these two to zero as much as possible
#define PARTICLE_SYSTEM_BASE_REAL 0
#define PARTICLE_SYSTEM_BASE_INT 0

class ParticleSystem
    : public amrex::NeighborParticleContainer<PARTICLE_SYSTEM_BASE_REAL, PARTICLE_SYSTEM_BASE_INT>
{
    private:

        //- Input file
        const ParmParse&    m_input;

    protected:

        /* GETTERS */

        // Access input parameters
        const ParamParse& input();

        //- Get the particle position (x-component)
        amrex::RealVector& posx(); 

        //- Get the particle position (y-component)
        amrex::RealVector& posy(); 

        //- Get the particle position (z-component)
        amrex::RealVector& posz(); 

        //- Get the particle velocity (x-component)
        amrex::RealVector& velx();

        //- Get the particle velocity (y-component)
        amrex::RealVector& vely();

        //- Get the particle velocity (z-component)
        amrex::RealVector& velz();

        //- Get the particle angular velocity (x-component)
        amrex::RealVector& avelx();

        //- Get the particle angular velocity (y-component)
        amrex::RealVector& avely();

        //- Get the particle angular velocity (z-component)
        amrex::RealVector& avelz();

        //- Get the particle radius
        amrex::RealVector& radius();

        //- Get the particle mass
        amrex::RealVector& mass();

        //- Get the particle id
        amrex::IntVector& pId();

        //- Get the particle processor
        amrex::IntVector& 

    public:

        ~ParticleSystem();

        //- Main constructor
        ParticleSystem
        ( 
            const ParmParse&                   input,
            const amrex::Geometry              geom,
            const amrex::DistributionMapping   dmap,
            const amrex::BoxArray              ba,
            int                                numcells,
        );

        /**
         * TODO: Add constructors with particle initialization
         */


};

