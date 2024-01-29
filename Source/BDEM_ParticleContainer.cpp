#include <BDEM_ParticleContainer.H>
#include <BDEM_Collide_Utils.H>
#include <stl_tools/STLtools.H>
#include <MoveUtils.H>
#include <WallTemp.H>

using namespace amrex;
using namespace std;

Real BDEMParticleContainer::compute_coll_timescale(int bonded_sphere_particles)
{
    BL_PROFILE("BDEMParticleContainer::compute_coll_timescale");
    Real coll_timescale=BIGVAL;

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);

    using PType = typename BDEMParticleContainer::SuperParticleType;
    coll_timescale = amrex::ReduceMin(*this, [=] 
    AMREX_GPU_HOST_DEVICE (const PType& p) -> Real 
    {
        //technically the mass should be effective mass 
        //among collision partners
        //but this may be ok for now
        // FIXME: should this change when using glued sphere model?
        // return( std::pow(DEM::k_n/p.rdata(realData::mass)
        Real bonded_dt = 1e10;
        if(bonded_sphere_particles) bonded_dt = 0.8165*2.0*p.rdata(realData::radius)*pow(p.rdata(realData::density) / DEM::E_bond,0.5); 
        return( std::min(bonded_dt, std::pow(DEM::k_n/(p.rdata(realData::mass))
                *PI*PI/(PI*PI + log(DEM::e_n)*log(DEM::e_n))
                ,-half))
         ); 
    });

#ifdef BL_USE_MPI
    ParallelDescriptor::ReduceRealMin(coll_timescale);
#endif

    return(coll_timescale);
}

void BDEMParticleContainer::computeForces (Real &dt,const EBFArrayBoxFactory *eb_factory, 
            const MultiFab *lsmfab,
            bool do_heat_transfer, int walltemp_vardir,
            Real walltemp_polynomial[3],
            const int ls_refinement,bool stl_geom_present, int contact_law, int steps,
            RealVect &gravity,
            const int bonded_sphere_particles,
            const int liquid_bridging,
            const int particle_cohesion, 
            const Real init_force, const int init_force_dir, const int init_force_comp,
            const Real cb_force, const Real cb_torq, const int cb_dir)
{
    BL_PROFILE("BDEMParticleContainer::computeForces");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    bool stlpresent=stl_geom_present;

    bool resolve_levset_wall_collisions=(eb_factory != NULL);
    
    std::map<PairIndex, bool> particle_tile_has_walls;

    if(resolve_levset_wall_collisions)
    {
        const FabArray<EBCellFlagFab>* flags = &(eb_factory->getMultiEBCellFlagFab());

        for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            int gid = mfi.index();
            int tid = mfi.LocalTileIndex();
            auto index = std::make_pair(gid, tid);
            const Box& bx = mfi.tilebox();

            Box phibx = convert(bx, {0,0,0});
            phibx.refine(ls_refinement);
            phibx.surroundingNodes();

            bool has_wall = false;
            if ((*flags)[mfi].getType(amrex::grow(bx,1)) == FabType::singlevalued)
            {
                particle_tile_has_walls[index] = true;
            }
        }
    }

    //zero forces
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        auto nbor_data = m_neighbor_list[lev][index].data();
        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
                AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            p.rdata(realData::fx) = zero;
            p.rdata(realData::fy) = zero;
            p.rdata(realData::fz) = zero;

            p.rdata(realData::taux) = zero;
            p.rdata(realData::tauy) = zero;
            p.rdata(realData::tauz) = zero;
        });
    }

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        auto nbor_data = m_neighbor_list[lev][index].data();
        ParticleType* pstruct = aos().dataPtr();
        
        int x_lo_bc = domain_bc[0];
        int x_hi_bc = domain_bc[1];
        int y_lo_bc = domain_bc[2];
        int y_hi_bc = domain_bc[3];
        int z_lo_bc = domain_bc[4];
        int z_hi_bc = domain_bc[5];

        if(resolve_levset_wall_collisions)
        {
            if (particle_tile_has_walls[index])
            {
#include"BDEM_LevsetWallCollisions.H"
            }
        }
        if(stlpresent)
        {
#include"BDEM_STLCollisions.H"
        }
        //at Cartesian walls
#include"BDEM_CartWallCollisions.H"
  

#include"BDEM_ParticleCollisions.H"
    }
}

void BDEMParticleContainer::computeForcing(Real time, Real time_offset, int forcing_type, Real forcing_freq,
                                          Real forcing_amp, Real forcing_dir[THREEDIM])
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        auto nbor_data = m_neighbor_list[lev][index].data();
        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
                  AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];

            Real forcing_f;
            if(forcing_type == 1){
                forcing_f = forcing_amp*sin(2*PI*(time+time_offset)*forcing_freq);
                p.rdata(realData::fx) += forcing_f*forcing_dir[0]*p.rdata(realData::mass);
                p.rdata(realData::fy) += forcing_f*forcing_dir[1]*p.rdata(realData::mass);
                p.rdata(realData::fz) += forcing_f*forcing_dir[2]*p.rdata(realData::mass);
            }
        });
    }

}

void BDEMParticleContainer::moveParticles(const amrex::Real& dt,
        int do_chemistry,Real minradfrac, int verlet_scheme)
{
    BL_PROFILE("BDEMParticleContainer::moveParticles");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    const auto dx = Geom(lev).CellSizeArray();
    auto& plev  = GetParticles(lev);

    //capture react stuff
    int nspecies=m_chemptr->nspecies;
    int nreac=m_chemptr->nreac;
    int nsolidspecs=m_chemptr->solidspec_ids.size();
    Real *reactmatrix=m_chemptr->reactmatrix.data();
    Real *arrh_A=m_chemptr->arrh_A.data();
    Real *arrh_Ea=m_chemptr->arrh_Ea.data();
    Real *molwts=m_chemptr->molwts.data();
    int *solidspec_ids=m_chemptr->solidspec_ids.data();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int x_lo_bc = domain_bc[0];
        int x_hi_bc = domain_bc[1];
        int y_lo_bc = domain_bc[2];
        int y_hi_bc = domain_bc[3];
        int z_lo_bc = domain_bc[4];
        int z_hi_bc = domain_bc[5];
        
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();
        ParticleType* pstruct = aos().dataPtr();

        // now we move the particles
        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            if(p.idata(intData::phase) != -1){    // Particles with phase = -1 are inert and do not move
              Real pos_old[3];
              Real pos_new[3];
              Real rp = p.rdata(realData::radius);

              pos_old[0]=p.pos(0);
              pos_old[1]=p.pos(1);
              pos_old[2]=p.pos(2);

              Real verlet_factor = (verlet_scheme) ? 0.5:1.0;

              // Global damping of particle forces
              p.rdata(realData::fx) -= DEM::global_damping * p.rdata(realData::xvel);
              p.rdata(realData::fy) -= DEM::global_damping * p.rdata(realData::yvel);
              p.rdata(realData::fz) -= DEM::global_damping * p.rdata(realData::zvel);

              // Force-based damping
              Real vel_vect[THREEDIM] = {p.rdata(realData::xvel), p.rdata(realData::yvel), p.rdata(realData::zvel)};
              Real f_vect[THREEDIM] = {p.rdata(realData::fx), p.rdata(realData::fy), p.rdata(realData::fz)};
              Real vmag = sqrt(dotpdt(vel_vect, vel_vect));
              Real fmag = sqrt(dotpdt(f_vect, f_vect));
              if(vmag > TINYVAL){
                  p.rdata(realData::fx) -= DEM::force_damping * fmag * (p.rdata(realData::xvel)/vmag);
                  p.rdata(realData::fy) -= DEM::force_damping * fmag * (p.rdata(realData::yvel)/vmag);
                  p.rdata(realData::fz) -= DEM::force_damping * fmag * (p.rdata(realData::zvel)/vmag);
              }

              p.rdata(realData::xvel) += (p.rdata(realData::fx)/p.rdata(realData::mass)) * dt * verlet_factor;
              p.rdata(realData::yvel) += (p.rdata(realData::fy)/p.rdata(realData::mass)) * dt * verlet_factor;
              p.rdata(realData::zvel) += (p.rdata(realData::fz)/p.rdata(realData::mass)) * dt * verlet_factor;

              if(verlet_scheme != 2){
                  p.pos(0) += p.rdata(realData::xvel) * dt;
                  p.pos(1) += p.rdata(realData::yvel) * dt;
                  p.pos(2) += p.rdata(realData::zvel) * dt;
              }

              // Global damping of particle torques
              p.rdata(realData::taux) -= pow(p.rdata(realData::radius),two)*DEM::global_damping*p.rdata(realData::xangvel);
              p.rdata(realData::tauy) -= pow(p.rdata(realData::radius),two)*DEM::global_damping*p.rdata(realData::yangvel);
              p.rdata(realData::tauz) -= pow(p.rdata(realData::radius),two)*DEM::global_damping*p.rdata(realData::zangvel);

              // Torque-based damping
              Real angvel_vect[THREEDIM] = {p.rdata(realData::xangvel), p.rdata(realData::yangvel), p.rdata(realData::zangvel)};
              Real tau_vect[THREEDIM] = {p.rdata(realData::taux), p.rdata(realData::tauy), p.rdata(realData::tauz)};
              Real angvmag = sqrt(dotpdt(angvel_vect, angvel_vect));
              Real tmag = sqrt(dotpdt(tau_vect, tau_vect));
              if(angvmag > TINYVAL){
                  p.rdata(realData::taux) -= DEM::force_damping * tmag * (p.rdata(realData::xangvel)/angvmag);
                  p.rdata(realData::tauy) -= DEM::force_damping * tmag * (p.rdata(realData::yangvel)/angvmag);
                  p.rdata(realData::tauz) -= DEM::force_damping * tmag * (p.rdata(realData::zangvel)/angvmag);
              }

              p.rdata(realData::xangvel) += p.rdata(realData::taux) * p.rdata(realData::Iinv) * dt * verlet_factor - DEM::angv_damping*p.rdata(realData::xangvel);
              p.rdata(realData::yangvel) += p.rdata(realData::tauy) * p.rdata(realData::Iinv) * dt * verlet_factor - DEM::angv_damping*p.rdata(realData::yangvel);
              p.rdata(realData::zangvel) += p.rdata(realData::tauz) * p.rdata(realData::Iinv) * dt * verlet_factor - DEM::angv_damping*p.rdata(realData::zangvel);

              // Tracking change in theta_x for beam twisting testing
              if(verlet_scheme != 2) p.rdata(realData::theta_x) += p.rdata(realData::xangvel) * dt;

              // FIXME: Chemistry should be compatible w Verlet scheme
              if(do_chemistry)
              {
                  Real wdot[MAXSPECIES+1]={0.0};
                  Real spec[MAXSPECIES]={0.0};
                  Real minrad=minradfrac*p.rdata(realData::radinit);

                  for(int sp=0;sp<nspecies;sp++)
                  {
                      //conc in kg/m3
                      spec[sp]=p.rdata(realData::firstspec+sp)*p.rdata(realData::density); 
                  }

                  getProductionRate(nspecies,nsolidspecs,nreac,spec,molwts,p.rdata(realData::density), 
                                    p.rdata(realData::radius), p.rdata(realData::radinit), p.rdata(realData::temperature),
                                    solidspec_ids, reactmatrix, arrh_A, arrh_Ea, wdot);

                  for(int sp=0;sp<nspecies;sp++)
                  {
                      p.rdata(realData::firstspec+sp) += wdot[sp]*dt/p.rdata(realData::density);
                  }
                  p.rdata(realData::radius) += wdot[nspecies]*dt;

                  //reset radius
                  if(p.rdata(realData::radius)<minrad)
                  {
                      p.rdata(realData::radius)=minrad;
                  }
              }

              // FIXME: Update for glued sphere code
              if (x_lo_bc==HARDWALL_BC and p.pos(0) < (plo[0]+rp))
              {
                  p.pos(0) = two*(plo[0]+rp) - p.pos(0);
                  p.rdata(realData::xvel) = -p.rdata(realData::xvel);
              }
              if (x_hi_bc==HARDWALL_BC and p.pos(0) > (phi[0]-rp))
              {
                  p.pos(0) = two*(phi[0]-rp) - p.pos(0);
                  p.rdata(realData::xvel) = -p.rdata(realData::xvel);
              }
              if (y_lo_bc==HARDWALL_BC and p.pos(1) < (plo[1]+rp))
              {
                  p.pos(1) = two*(plo[1]+rp) - p.pos(1);
                  p.rdata(realData::yvel) = -p.rdata(realData::yvel);
              }
              if (y_hi_bc==HARDWALL_BC and p.pos(1) > (phi[1]-rp))
              {
                  p.pos(1) = two*(phi[1]-rp) - p.pos(1);
                  p.rdata(realData::yvel) = -p.rdata(realData::yvel);
              }
              if (z_lo_bc==HARDWALL_BC and p.pos(2) < (plo[2]+rp))
              {
                  p.pos(2) = two*(plo[2]+rp) - p.pos(2);
                  p.rdata(realData::zvel) = -p.rdata(realData::zvel);
              }
              if (z_hi_bc==HARDWALL_BC and p.pos(2) > (phi[2]-rp))
              {
                  p.pos(2) = two*(phi[2]-rp) - p.pos(2);
                  p.rdata(realData::zvel) = -p.rdata(realData::zvel);
              }
            }

        });

    }
}

void BDEMParticleContainer::clipParticles(int clip_particle_dir, Real clip_particle_val)
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        ParticleType* pstruct = aos().dataPtr();
        const Box & bx = mfi.tilebox();
        amrex::ParallelFor(np,[=] AMREX_GPU_DEVICE (int i) noexcept{
            ParticleType& p = pstruct[i];
            Real ppos_inert[THREEDIM];
            ppos_inert[XDIR] = p.pos(0);
            ppos_inert[YDIR] = p.pos(1);
            ppos_inert[ZDIR] = p.pos(2);

            if(ppos_inert[clip_particle_dir] > clip_particle_val)
            {
                p.id()=-1;
            }
        });
    }
    Redistribute();
}

void BDEMParticleContainer::saveParticles_softwall()
{
    BL_PROFILE("BDEMParticleContainer::saveParticles_softwall");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    auto& plev  = GetParticles(lev);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int x_lo_bc = domain_bc[0];
        int x_hi_bc = domain_bc[1];
        int y_lo_bc = domain_bc[2];
        int y_hi_bc = domain_bc[3];
        int z_lo_bc = domain_bc[4];
        int z_hi_bc = domain_bc[5];
        
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();
        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            p.rdata(realData::posx_prvs)=p.pos(0);
            p.rdata(realData::posy_prvs)=p.pos(1);
            p.rdata(realData::posz_prvs)=p.pos(2);
            p.idata(intData::near_softwall)=0;
            
            if (x_lo_bc==SOFTWALL_BC and p.pos(0) < plo[0])
            {
                p.pos(0) = two*plo[0] - p.pos(0);
                p.idata(intData::near_softwall)=1;
            }
            if (x_hi_bc==SOFTWALL_BC and p.pos(0) > phi[0])
            {
                p.pos(0) = two*phi[0] - p.pos(0);
                p.idata(intData::near_softwall)=1;
            }
            if (y_lo_bc==SOFTWALL_BC and p.pos(1) < plo[1])
            {
                p.pos(1) = two*plo[1] - p.pos(1);
                p.idata(intData::near_softwall)=1;
            }
            if (y_hi_bc==SOFTWALL_BC and p.pos(1) > phi[1])
            {
                p.pos(1) = two*phi[1] - p.pos(1);
                p.idata(intData::near_softwall)=1;
            }
            if (z_lo_bc==SOFTWALL_BC and p.pos(2) < plo[2])
            {
                p.pos(2) = two*plo[2] - p.pos(2);
                p.idata(intData::near_softwall)=1;
            }
            if (z_hi_bc==SOFTWALL_BC and p.pos(2) > phi[2])
            {
                p.pos(2) = two*phi[2] - p.pos(2);
                p.idata(intData::near_softwall)=1;
            }
        });
    }
}

void BDEMParticleContainer::reassignParticles_softwall()
{
    BL_PROFILE("BDEMParticleContainer::reassignParticles_softwall");

    const int lev = 0;
    auto& plev  = GetParticles(lev);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();
        ParticleType* pstruct = aos().dataPtr();

        // now we move the particles
        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            if(p.idata(intData::near_softwall)==1)
            {
                p.pos(0)=p.rdata(realData::posx_prvs);
                p.pos(1)=p.rdata(realData::posy_prvs);
                p.pos(2)=p.rdata(realData::posz_prvs);
                p.idata(intData::near_softwall)=0;
            }
        });
    }
}

void BDEMParticleContainer::computeMoistureContent(Real MC_avg, Real MC_stdev, Real liquid_density, Real FSP)
{
    BL_PROFILE("BDEMParticleContainer::computeMoistureContent");

    const int lev = 0;
    auto& plev = GetParticles(lev);

    // Calculate the liquid volume per particle (assumed proportional to diameter^2) 
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();
        ParticleType* pstruct = aos().dataPtr();

        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            // Real MC = (MC_stdev > 0.0) ? std::min(std::max(amrex::RandomNormal(MC_avg, MC_stdev),0.0),0.9):MC_avg;
            Real MC = MC_avg;
	    Real original_lv = p.rdata(realData::liquid_volume);
            p.rdata(realData::liquid_volume) = (MC > FSP || MC == 0.0) ? (p.rdata(realData::density)*p.rdata(realData::volume)/liquid_density)*(MC - FSP)/(1.0 - MC):0;
            p.rdata(realData::mass) = p.rdata(realData::density)*p.rdata(realData::volume) * (1.0 + MC/(1.0 - MC));
            p.rdata(realData::density) = p.rdata(realData::mass) / p.rdata(realData::volume);

	    // Rescale the total liquid participating in bridges
	    if(original_lv > 0.0) p.rdata(realData::total_bridge_volume) *= p.rdata(realData::liquid_volume) / original_lv;

	    // If MC < FSP, remove all existing liquid bridges
	    for(int br=0; br<MAXBRIDGES; br++) p.idata(intData::first_bridge+br) = -1;
        });
    }
}

void BDEMParticleContainer::Calculate_Total_Mass_MaterialPoints(Real &total_mass, int cdir, Real cutoff)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    total_mass=0.0;

    using PType = typename BDEMParticleContainer::SuperParticleType;
    total_mass = amrex::ReduceSum(*this, [=]
        AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real ppos[3] = {p.pos(0), p.pos(1), p.pos(2)};
            if(p.idata(intData::phase)==0 && ppos[cdir] > cutoff)
            {
              return(p.rdata(realData::mass));
            }
            else
            {
              return(0.0);
            }
        });

  #ifdef BL_USE_MPI
      ParallelDescriptor::ReduceRealSum(total_mass);
  #endif
}

void BDEMParticleContainer::Calculate_Total_Speed_MaterialPoints(Real &total_speed)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto domain = geom.Domain();

    total_speed=0.0;

    using PType = typename BDEMParticleContainer::SuperParticleType;
    total_speed = amrex::ReduceSum(*this, [=]
        AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            Real pvel[3] = {p.rdata(realData::xvel), p.rdata(realData::yvel), p.rdata(realData::zvel)};
            return(sqrt(pvel[0]*pvel[0] + pvel[1]*pvel[1] + pvel[2]*pvel[2]));
        });

  #ifdef BL_USE_MPI
      ParallelDescriptor::ReduceRealSum(total_speed);
  #endif
}

void BDEMParticleContainer::writeParticles(const int n, const int bonded_sphere_particles, const std::string pltprefix)
{
    BL_PROFILE("BDEMParticleContainer::writeParticles");
    const std::string& pltfile = pltprefix + amrex::Concatenate("plt", n, 5);

    Vector<int> writeflags_real(realData::count+MAXSPECIES-1,0);
    Vector<int> writeflags_int(intData::count,0);

    Vector<std::string> real_data_names;
    Vector<std::string>  int_data_names;

    real_data_names.push_back("radius");
    real_data_names.push_back("radinit");
    real_data_names.push_back("xvel");
    real_data_names.push_back("yvel");
    real_data_names.push_back("zvel");
    real_data_names.push_back("fx");
    real_data_names.push_back("fy");
    real_data_names.push_back("fz");
    real_data_names.push_back("xangvel");
    real_data_names.push_back("yangvel");
    real_data_names.push_back("zangvel");
    real_data_names.push_back("taux");
    real_data_names.push_back("tauy");
    real_data_names.push_back("tauz");
    real_data_names.push_back("Iinv");
    real_data_names.push_back("volume");
    real_data_names.push_back("mass");
    real_data_names.push_back("density");
    real_data_names.push_back("E");
    real_data_names.push_back("nu");
    real_data_names.push_back("temperature");
    real_data_names.push_back("posx_prvs");
    real_data_names.push_back("posy_prvs");
    real_data_names.push_back("posz_prvs");
    real_data_names.push_back("euler_angle_x");
    real_data_names.push_back("euler_angle_y");
    real_data_names.push_back("euler_angle_z");
    real_data_names.push_back("liquid_volume");
    real_data_names.push_back("total_bridge_volume");
    real_data_names.push_back("theta_x");

    // For debug
    // real_data_names.push_back("overlap_n");
    // real_data_names.push_back("kn");
    // real_data_names.push_back("kt");
    // real_data_names.push_back("eta_n");
    // real_data_names.push_back("eta_t");
    // real_data_names.push_back("overlap_n_y");
    // real_data_names.push_back("fn_y");
    // real_data_names.push_back("vr_n_y");
    // real_data_names.push_back("overlap_t_y");
    // real_data_names.push_back("ft_y");
    // real_data_names.push_back("vr_t_y");

    for(int i=0;i<MAXBONDS;i++)
    {
        std::string bondval = amrex::Concatenate("bval_fx",i,2);
        real_data_names.push_back(bondval);
        bondval = amrex::Concatenate("bval_fy",i,2);
        real_data_names.push_back(bondval);
        bondval = amrex::Concatenate("bval_fz",i,2);
        real_data_names.push_back(bondval);
        bondval = amrex::Concatenate("bval_tnx",i,2);
        real_data_names.push_back(bondval);
        bondval = amrex::Concatenate("bval_tny",i,2);
        real_data_names.push_back(bondval);
        bondval = amrex::Concatenate("bval_tnz",i,2);
        real_data_names.push_back(bondval);
        bondval = amrex::Concatenate("bval_ttx",i,2);
        real_data_names.push_back(bondval);
        bondval = amrex::Concatenate("bval_tty",i,2);
        real_data_names.push_back(bondval);
        bondval = amrex::Concatenate("bval_ttz",i,2);
        real_data_names.push_back(bondval);
    }

    for(int i=0;i<MAXSPECIES;i++)
    {
        std::string chemname = amrex::Concatenate("chemspec",i,2);
        real_data_names.push_back(chemname);
    }

    for(int i=0;i<m_chemptr->nspecies;i++)
    {
        real_data_names[realData::firstspec+i]=m_chemptr->specnames[i];
    }

    int_data_names.push_back("phase");
    int_data_names.push_back("near_softwall");
    int_data_names.push_back("type_id");

    for(int i=0;i<MAXBRIDGES;i++)
    {
        std::string bridgeidx = amrex::Concatenate("p2idx",i,2);
        int_data_names.push_back(bridgeidx);
    }

    for(int i=0;i<MAXBONDS;i++)
    {
        std::string bondidx = amrex::Concatenate("bidx",i,2);
        int_data_names.push_back(bondidx);
    }

    writeflags_real[realData::radius]=1;
    writeflags_real[realData::xvel]=1;
    writeflags_real[realData::yvel]=1;
    writeflags_real[realData::zvel]=1;
    writeflags_real[realData::fx]=1;
    writeflags_real[realData::fy]=1;
    writeflags_real[realData::fz]=1;
    writeflags_real[realData::xangvel]=1;
    writeflags_real[realData::yangvel]=1;
    writeflags_real[realData::zangvel]=1;
    writeflags_real[realData::taux]=1;
    writeflags_real[realData::tauy]=1;
    writeflags_real[realData::tauz]=1;
    writeflags_real[realData::mass]=1;
    writeflags_real[realData::temperature]=1;
    writeflags_real[realData::liquid_volume]=1;
    writeflags_real[realData::total_bridge_volume]=1;
    writeflags_real[realData::theta_x]=1;

    for(int i=0;i<m_chemptr->nspecies;i++)
    {
        writeflags_real[realData::firstspec+i]=1;
    }
    if(bonded_sphere_particles){
        writeflags_int[intData::type_id] = 1;
    }

    WritePlotFile(pltfile, "particles",writeflags_real, 
            writeflags_int, real_data_names, int_data_names);

}

void BDEMParticleContainer::removeEBOverlapParticles(EBFArrayBoxFactory *eb_factory,
                                                     const MultiFab *lsmfab, const int ls_refinement){
  // This function removes particles that have gotten "wedged" in corners of an EB
  // boundary (where forces are not calculated correctly, and particles may stick out of the boundary)

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();

    bool resolve_levset_wall_collisions=(eb_factory != NULL);
    
    std::map<PairIndex, bool> particle_tile_has_walls;

    if(resolve_levset_wall_collisions)
    {
        const FabArray<EBCellFlagFab>* flags = &(eb_factory->getMultiEBCellFlagFab());

        for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            int gid = mfi.index();
            int tid = mfi.LocalTileIndex();
            auto index = std::make_pair(gid, tid);
            const Box& bx = mfi.tilebox();

            Box phibx = convert(bx, {0,0,0});
            phibx.refine(ls_refinement);
            phibx.surroundingNodes();

            bool has_wall = false;
            if ((*flags)[mfi].getType(amrex::grow(bx,1)) == FabType::singlevalued)
            {
                auto& ptile = plev[index];
                auto& aos   = ptile.GetArrayOfStructs();
                const size_t np = aos.numParticles();

                ParticleType* pstruct = aos().dataPtr();

                const auto phiarr = lsmfab->array(mfi);

                amrex::ParallelFor(np,[=] AMREX_GPU_DEVICE (int i) noexcept{
                    ParticleType& p = pstruct[i];
                    Real ppos_inert[THREEDIM];
                    ppos_inert[XDIR] = p.pos(0);
                    ppos_inert[YDIR] = p.pos(1);
                    ppos_inert[ZDIR] = p.pos(2);

                    Real ls_value = get_levelset_value(ppos_inert, ls_refinement, phiarr, plo, dx);

                    if (ls_value < p.rdata(realData::radius)*0.95)
                    {
                        p.id() = -1;
                    }
                });
            }
        }
    }
    Redistribute();
}



