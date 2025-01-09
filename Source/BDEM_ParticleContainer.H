#ifndef BDEM_PARTICLE_CONTAINER_H_
#define BDEM_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>
#include <BDEM_Specs.H>
#include <Constants.H>
#include <Mechanism.H>
#include <BDEM_CheckPair.H>
#include <BDEM_ParticleTypes.H>
#include <stl_tools/STLtools.H>


class BDEMParticleContainer
    : public amrex::NeighborParticleContainer<realData::count+MAXSPECIES-1, intData::count>
{

public:

    BDEMParticleContainer (const amrex::Geometry            geom,
                         const amrex::DistributionMapping   dmap,
                         const amrex::BoxArray              ba,
                         int                                numcells,
                         mechanism *chemptr)

    : NeighborParticleContainer<realData::count+MAXSPECIES-1, intData::count>(geom, dmap, ba, numcells)
    {
        m_chemptr=chemptr;

        // Particle properties that need to be copied for neighbor updates
        setRealCommComp(0, true);   // posx
        setRealCommComp(1, true);   // posy
        setRealCommComp(2, true);   // posz
        setRealCommComp(3+realData::mass, true); 
        setRealCommComp(3+realData::radius, true); 
        setRealCommComp(3+realData::xvel, true); 
        setRealCommComp(3+realData::yvel, true); 
        setRealCommComp(3+realData::zvel, true); 
        setRealCommComp(3+realData::xangvel, true); 
        setRealCommComp(3+realData::yangvel, true); 
        setRealCommComp(3+realData::zangvel, true); 
        setRealCommComp(3+realData::liquid_volume, true); 
        setRealCommComp(3+realData::E, true); 
        setRealCommComp(3+realData::nu, true); 
        setRealCommComp(3+realData::temperature, true); 

        // // Remaining particle properties don't need to be copied over
        setRealCommComp(3+realData::radinit, false); 
        setRealCommComp(3+realData::fx, false); 
        setRealCommComp(3+realData::fy, false); 
        setRealCommComp(3+realData::fz, false); 
        setRealCommComp(3+realData::taux, false);
        setRealCommComp(3+realData::tauy, false);
        setRealCommComp(3+realData::tauz, false);
        setRealCommComp(3+realData::Iinv, false);
        setRealCommComp(3+realData::volume, false);
        setRealCommComp(3+realData::density, false);
        setRealCommComp(3+realData::posx_prvs, false);
        setRealCommComp(3+realData::posy_prvs, false);
        setRealCommComp(3+realData::posz_prvs, false);
        setRealCommComp(3+realData::euler_angle_x, false);
        setRealCommComp(3+realData::euler_angle_y, false);
        setRealCommComp(3+realData::euler_angle_z, false);
        setRealCommComp(3+realData::total_bridge_volume, false);
        for(int i=0; i<MAXBONDS*9; i++) setRealCommComp(3+realData::first_bond_v+i, false);
        setRealCommComp(3+realData::theta_x, false);
        for(int i=0; i<MAXSPECIES; i++) setRealCommComp(3+realData::firstspec+i, false);
        setRealCommComp(3+realData::fraction_of_fibrils, false); 
        // setIntCommComp(intData::phase, false);
        // setIntCommComp(intData::near_softwall, false);
        // setIntCommComp(intData::type_id, false);
        // for(int i=0; i<MAXBRIDGES; i++) setIntCommComp(intData::first_bridge+i, false);
        // for(int i=0; i<MAXBONDS; i++) setIntCommComp(intData::first_bond+i, false);
    }

    void InitParticles (const std::string & filename,
                        bool &do_heat_transfer, 
                        Real temp_mean, Real temp_stdev,
                        int contact_law, int liquid_bridging,
                        Real liquid_density, Real MC_avg, 
                        Real MC_stdev, Real FSP, const int solve_fibrils = 0);

    void InitParticles (Real mincoords[THREEDIM],Real maxcoords[THREEDIM], 
        Real meanvel[THREEDIM], Real fluctation[THREEDIM], Real rad, 
        Real dens, Real E, Real nu, Real temp, Real spec[MAXSPECIES], 
        int do_multi_part_per_cell, int layer_particles,
        int bonded_sphere_particles,Real min_rad, Real max_rad,  
        int p_type, Vector<int> type_list, int use_type_dist, 
        Vector<Real> dist_list,int liquid_bridging, 
        Real liquid_density, Real MC_avg, Real MC_stdev, Real FSP
        );

    void clearSourcingVolume(Real mincoords[THREEDIM],Real maxcoords[THREEDIM]);

    void InitBondedParticles (const std::string & filename,
                              bool &do_heat_transfer,
                              int contact_law,
                              int cantilever_beam_test,
                              int liquid_bridging,
                              Real liquid_density, Real MC_avg, 
                              Real MC_stdev, Real FSP);

    void InitChemSpecies(amrex::Real Yis[MAXSPECIES]);

    void InitChemSpecies(int ndomains, Real *mincoords,
        Real *maxcoords,Real *spec_massfracs);

    void computeForces (Real &dt,const EBFArrayBoxFactory *eb_factory, 
            const MultiFab *lsmfab,
            bool do_heat_transfer, int walltemp_vardir,
            Real walltemp_polynomial[3],
            const int ls_refinement,bool stl_geom_resent, int contact_law, int steps,
            RealVect &gravity, amrex::Vector< stl_specs >& stls, Real time,
            const int bonded_sphere_particles=0,
            const int liquid_bridging=0, 
            const int particle_cohesion=0,
            const Real init_force = 0.0, const int init_force_dir = 0, const int init_force_comp = 0,
            const Real cb_force = 0.0, const Real cb_torq = 0.0, const int cb_dir = 0, const int drag_model = 0,
            const int solve_fibrillation = 0);

    Real compute_coll_timescale(int bonded_sphere_particles);

    void removeParticlesOutsideBoundary(const MultiFab *lsfab,const EBFArrayBoxFactory *ebfactory,
            const int ls_refinement);
    
    void removeParticlesInsideSTL(amrex::Vector<amrex::Real> outside_point, STLtools* stl);
    void checkParticlesInsideSTL(amrex::Vector<amrex::Real> outside_point, STLtools* stl);

    void moveParticles (const amrex::Real& dt, int do_chemistry=0, amrex::Real minradfrac=0.1, int verlet_schme = 0);

    void writeParticles (const int n, const int bonded_sphere_particles, const std::string pltprefix);   

    ParticleType generate_particle(Real x,Real y,Real z,
                Real velx, Real vely, Real velz,
                Real dens, Real E, Real nu, Real rad, Real temp,Real spec[MAXSPECIES], 
                Real liquid_density, Real MC, Real FSP, int p_type=0);

    void redist_particles(const int local,bool using_softwalls)
    {
        const int lev_min = 0;
        const int lev_max = 0;
        const int nGrow = 0;
        if(using_softwalls)
        {
            saveParticles_softwall(); 
        }
        Redistribute (lev_min,lev_max,nGrow,local);

        fillNeighbors();
        buildNeighborList(CheckPair());
    }

    void set_domain_bcs(Vector<int> bclo,Vector<int> bchi)
    {
        domain_bc[0]=bclo[0];
        domain_bc[1]=bchi[0];

        domain_bc[2]=bclo[1];
        domain_bc[3]=bchi[1];

        domain_bc[4]=bclo[2];
        domain_bc[5]=bchi[2];
    }


    void update_reflectscore();
    void reflect_off_moving_stl();

    void saveParticles_softwall();
    void reassignParticles_softwall(); 

    void computeMoistureContent(Real MC_avg, Real MC_stdev, Real liquid_density, Real FSP);

    void clipParticles(int clip_particle_dir, Real clip_particle_val);
    
    void removeEBOverlapParticles(EBFArrayBoxFactory *eb_factory, const MultiFab *lsmfab, const int ls_refinement);

    void Calculate_Total_Mass_MaterialPoints(Real &total_mass, int cdir, Real cutoff);

    void Calculate_Total_Speed_MaterialPoints(Real &total_speed);

    void reassignParticleProperties(Real reinit_rad, Real reinit_dens, Real reinit_E, Real reinit_nu);

    void particleColoringStriation(Real striation_len, int striation_dir);

private:

    int domain_bc[6]={0};
    mechanism *m_chemptr;
};
#endif
