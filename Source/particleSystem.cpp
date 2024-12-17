#include "particleSystem.h"


/* CONSTRUCTORS AND DESRUCTOR */
ParticleSystem::ParticleSystem
( 
    const ParmParse&                   input,
    const amrex::Geometry              geom,
    const amrex::DistributionMapping   dmap,
    const amrex::BoxArray              ba,
    int                                numcells,
)
:
    amrex::NeighborParticleContainer<PARTICLE_SYSTEM_BASE_REAL, PARTICLE_SYSTEM_BASE_INT>(input, geom, dmap, ba, numcells),
    m_input(input)
{

    //- Add linear velocity, angular velocity, radius, and mass
    AddRealComp("velx",true);
    AddRealComp("vely",true);
    AddRealComp("velz",true);
    AddRealComp("avelx",true);
    AddRealComp("avely",true);
    AddRealComp("avelz",true);
    AddRealComp("radius",true);
    AddRealComp("mass",true);
}

ParticleSystem::~ParticleSystem()
{
}

/* GETTERS */
const ParmParse& ParticleSystem::input()
{
    return m_input;
}

amrex::RealVector&  ParticleSystem::posx()
{
    return amrex::getRealData()[0];
}

amrex::RealVector&  ParticleSystem::posy()
{
    return amrex::getRealData()[1];
}

amrex::RealVector&  ParticleSystem::posz()
{
    return amrex::getRealData()[2];
}

amrex::RealVector&  ParticleSystem::velx()
{
    return amrex::getRealData("velx");
}

amrex::RealVector&  ParticleSystem::vely()
{
    return amrex::getRealData("vely");
}

amrex::RealVector&  ParticleSystem::velz()
{
    return amrex::getRealData("velz");
}

amrex::RealVector&  ParticleSystem::avelx()
{
    return amrex::getRealData("avelx");
}

amrex::RealVector&  ParticleSystem::avely()
{
    return amrex::getRealData("avely");
}

amrex::RealVector&  ParticleSystem::avelz()
{
    return amrex::getRealData("avelz");
}

amrex::RealVector&  ParticleSystem::radius()
{
    return amrex::getRealData("radius");
}

amrex::RealVector&  ParticleSystem::mass()
{
    return amrex::getRealData("mass");
}
