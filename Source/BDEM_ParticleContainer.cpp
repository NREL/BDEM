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
        if(bonded_sphere_particles) bonded_dt = 0.8165*2.0*p.rdata(aos_realData::radius)*pow(p.rdata(aos_realData::density) / DEM::E_bond,0.5); 
        return( std::min(bonded_dt, std::pow(DEM::k_n/(p.rdata(aos_realData::mass))
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
        auto& soa   = ptile.GetStructOfArrays();
        const size_t np = aos.numParticles();

        auto nbor_data = m_neighbor_list[lev][index].data();
        ParticleType* pstruct = aos().dataPtr();

        auto fx_arr = soa.GetRealData(soa_realData::fx).data();
        auto fy_arr = soa.GetRealData(soa_realData::fy).data();
        auto fz_arr = soa.GetRealData(soa_realData::fz).data();

        auto taux_arr = soa.GetRealData(soa_realData::taux).data();
        auto tauy_arr = soa.GetRealData(soa_realData::tauy).data();
        auto tauz_arr = soa.GetRealData(soa_realData::tauz).data();

        amrex::ParallelFor(np,[=]
                AMREX_GPU_DEVICE (int i) noexcept
        {
            fx_arr[i] = zero;
            fy_arr[i] = zero;
            fz_arr[i] = zero;

            taux_arr[i] = zero;
            tauy_arr[i] = zero;
            tauz_arr[i] = zero;
        });
    }

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        auto& soa   = ptile.GetStructOfArrays();
        const size_t np = aos.numParticles();

        auto nbor_data = m_neighbor_list[lev][index].data();
        ParticleType* pstruct = aos().dataPtr();

        // Get SoA data for variables that are updated during force calculation kernels 
        auto fx_arr = soa.GetRealData(soa_realData::fx).data();
        auto fy_arr = soa.GetRealData(soa_realData::fy).data();
        auto fz_arr = soa.GetRealData(soa_realData::fz).data();

        auto taux_arr = soa.GetRealData(soa_realData::taux).data();
        auto tauy_arr = soa.GetRealData(soa_realData::tauy).data();
        auto tauz_arr = soa.GetRealData(soa_realData::tauz).data();

        auto total_bridge_vol_arr = soa.GetRealData(soa_realData::total_bridge_volume).data();
        Vector<int*> bridges_vec(MAXBRIDGES);
        for(int br=0; br<MAXBRIDGES; br++) bridges_vec[br] = soa.GetIntData(soa_intData::first_bridge+br).data(); 
        Vector<Real*> bonds_vec(MAXBONDS*9);
        for(int b=0; b<MAXBONDS*9; b++) bonds_vec[b] = soa.GetRealData(soa_realData::first_bond_v+b).data(); 
        
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
        auto& soa   = ptile.GetStructOfArrays();
        const size_t np = aos.numParticles();
        ParticleType* pstruct = aos().dataPtr();

        auto fx_arr = soa.GetRealData(soa_realData::fx).data();
        auto fy_arr = soa.GetRealData(soa_realData::fy).data();
        auto fz_arr = soa.GetRealData(soa_realData::fz).data();

        auto taux_arr = soa.GetRealData(soa_realData::taux).data();
        auto tauy_arr = soa.GetRealData(soa_realData::tauy).data();
        auto tauz_arr = soa.GetRealData(soa_realData::tauz).data();
        auto Iinv_arr = soa.GetRealData(soa_realData::Iinv).data();

        auto radinit_arr = soa.GetRealData(soa_realData::radinit).data();
        Vector<Real*> spec_vec(MAXSPECIES);
        for(int sp=0; sp<MAXSPECIES; sp++) spec_vec[sp] = soa.GetRealData(soa_realData::firstspec+sp).data(); 

        // now we move the particles
        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            if(soa.GetIntData(soa_intData::phase).data()[i] != -1){    // Particles with phase = -1 are inert and do not move
              Real pos_old[3];
              Real pos_new[3];
              Real rp = p.rdata(aos_realData::radius);

              pos_old[0]=p.pos(0);
              pos_old[1]=p.pos(1);
              pos_old[2]=p.pos(2);

              Real verlet_factor = (verlet_scheme) ? 0.5:1.0;

              // Global damping of particle forces
              fx_arr[i] -= DEM::global_damping * p.rdata(aos_realData::xvel);
              fy_arr[i] -= DEM::global_damping * p.rdata(aos_realData::yvel);
              fz_arr[i] -= DEM::global_damping * p.rdata(aos_realData::zvel);

              // Force-based damping
              Real vel_vect[THREEDIM] = {p.rdata(aos_realData::xvel), p.rdata(aos_realData::yvel), p.rdata(aos_realData::zvel)};
              Real f_vect[THREEDIM] = {fx_arr[i], fy_arr[i], fz_arr[i]};
              Real vmag = sqrt(dotpdt(vel_vect, vel_vect));
              Real fmag = sqrt(dotpdt(f_vect, f_vect));
              if(vmag > TINYVAL){
                  fx_arr[i] -= DEM::force_damping * fmag * (p.rdata(aos_realData::xvel)/vmag);
                  fy_arr[i] -= DEM::force_damping * fmag * (p.rdata(aos_realData::yvel)/vmag);
                  fz_arr[i] -= DEM::force_damping * fmag * (p.rdata(aos_realData::zvel)/vmag);
              }

              p.rdata(aos_realData::xvel) += (fx_arr[i]/p.rdata(aos_realData::mass)) * dt * verlet_factor;
              p.rdata(aos_realData::yvel) += (fy_arr[i]/p.rdata(aos_realData::mass)) * dt * verlet_factor;
              p.rdata(aos_realData::zvel) += (fz_arr[i]/p.rdata(aos_realData::mass)) * dt * verlet_factor;

              if(verlet_scheme != 2){
                  p.pos(0) += p.rdata(aos_realData::xvel) * dt;
                  p.pos(1) += p.rdata(aos_realData::yvel) * dt;
                  p.pos(2) += p.rdata(aos_realData::zvel) * dt;
              }

              // Global damping of particle torques
              taux_arr[i] -= pow(p.rdata(aos_realData::radius),two)*DEM::global_damping*p.rdata(aos_realData::xangvel);
              tauy_arr[i] -= pow(p.rdata(aos_realData::radius),two)*DEM::global_damping*p.rdata(aos_realData::yangvel);
              tauz_arr[i] -= pow(p.rdata(aos_realData::radius),two)*DEM::global_damping*p.rdata(aos_realData::zangvel);

              // Torque-based damping
              Real angvel_vect[THREEDIM] = {p.rdata(aos_realData::xangvel), p.rdata(aos_realData::yangvel), p.rdata(aos_realData::zangvel)};
              Real tau_vect[THREEDIM] = {taux_arr[i], tauy_arr[i], tauz_arr[i]};
              Real angvmag = sqrt(dotpdt(angvel_vect, angvel_vect));
              Real tmag = sqrt(dotpdt(tau_vect, tau_vect));
              if(angvmag > TINYVAL){
                  taux_arr[i] -= DEM::force_damping * tmag * (p.rdata(aos_realData::xangvel)/angvmag);
                  tauy_arr[i] -= DEM::force_damping * tmag * (p.rdata(aos_realData::yangvel)/angvmag);
                  tauz_arr[i] -= DEM::force_damping * tmag * (p.rdata(aos_realData::zangvel)/angvmag);
              }

              p.rdata(aos_realData::xangvel) += taux_arr[i] * Iinv_arr[i] * dt * verlet_factor - DEM::angv_damping*p.rdata(aos_realData::xangvel);
              p.rdata(aos_realData::yangvel) += tauy_arr[i] * Iinv_arr[i] * dt * verlet_factor - DEM::angv_damping*p.rdata(aos_realData::yangvel);
              p.rdata(aos_realData::zangvel) += tauz_arr[i] * Iinv_arr[i] * dt * verlet_factor - DEM::angv_damping*p.rdata(aos_realData::zangvel);

              // FIXME: Chemistry should be compatible w Verlet scheme
              if(do_chemistry)
              {
                  Real wdot[MAXSPECIES+1]={0.0};
                  Real spec[MAXSPECIES]={0.0};
                  Real minrad=minradfrac*radinit_arr[i];

                  for(int sp=0;sp<nspecies;sp++)
                  {
                      //conc in kg/m3
                      spec[sp]=spec_vec[sp][i]*p.rdata(aos_realData::density); 
                  }

                  getProductionRate(nspecies,nsolidspecs,nreac,spec,molwts,p.rdata(aos_realData::density), 
                                    p.rdata(aos_realData::radius), radinit_arr[i], p.rdata(aos_realData::temperature),
                                    solidspec_ids, reactmatrix, arrh_A, arrh_Ea, wdot);

                  for(int sp=0;sp<nspecies;sp++)
                  {
                      spec_vec[sp][i] += wdot[sp]*dt/p.rdata(aos_realData::density);
                  }
                  p.rdata(aos_realData::radius) += wdot[nspecies]*dt;

                  //reset radius
                  if(p.rdata(aos_realData::radius)<minrad)
                  {
                      p.rdata(aos_realData::radius)=minrad;
                  }
              }

              // FIXME: Update for glued sphere code
              if (x_lo_bc==HARDWALL_BC and p.pos(0) < (plo[0]+rp))
              {
                  p.pos(0) = two*(plo[0]+rp) - p.pos(0);
                  p.rdata(aos_realData::xvel) = -p.rdata(aos_realData::xvel);
              }
              if (x_hi_bc==HARDWALL_BC and p.pos(0) > (phi[0]-rp))
              {
                  p.pos(0) = two*(phi[0]-rp) - p.pos(0);
                  p.rdata(aos_realData::xvel) = -p.rdata(aos_realData::xvel);
              }
              if (y_lo_bc==HARDWALL_BC and p.pos(1) < (plo[1]+rp))
              {
                  p.pos(1) = two*(plo[1]+rp) - p.pos(1);
                  p.rdata(aos_realData::yvel) = -p.rdata(aos_realData::yvel);
              }
              if (y_hi_bc==HARDWALL_BC and p.pos(1) > (phi[1]-rp))
              {
                  p.pos(1) = two*(phi[1]-rp) - p.pos(1);
                  p.rdata(aos_realData::yvel) = -p.rdata(aos_realData::yvel);
              }
              if (z_lo_bc==HARDWALL_BC and p.pos(2) < (plo[2]+rp))
              {
                  p.pos(2) = two*(plo[2]+rp) - p.pos(2);
                  p.rdata(aos_realData::zvel) = -p.rdata(aos_realData::zvel);
              }
              if (z_hi_bc==HARDWALL_BC and p.pos(2) > (phi[2]-rp))
              {
                  p.pos(2) = two*(phi[2]-rp) - p.pos(2);
                  p.rdata(aos_realData::zvel) = -p.rdata(aos_realData::zvel);
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
        auto& soa   = ptile.GetStructOfArrays();
        const size_t np = aos.numParticles();
        ParticleType* pstruct = aos().dataPtr();

        auto pxp_arr = soa.GetRealData(soa_realData::posx_prvs).data();
        auto pyp_arr = soa.GetRealData(soa_realData::posy_prvs).data();
        auto pzp_arr = soa.GetRealData(soa_realData::posz_prvs).data();

        auto softwall_arr = soa.GetIntData(soa_intData::near_softwall).data();

        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            pxp_arr[i]=p.pos(0);
            pyp_arr[i]=p.pos(1);
            pzp_arr[i]=p.pos(2);
            softwall_arr[i]=0;
            
            if (x_lo_bc==SOFTWALL_BC and p.pos(0) < plo[0])
            {
                p.pos(0) = two*plo[0] - p.pos(0);
                softwall_arr[i]=1;
            }
            if (x_hi_bc==SOFTWALL_BC and p.pos(0) > phi[0])
            {
                p.pos(0) = two*phi[0] - p.pos(0);
                softwall_arr[i]=1;
            }
            if (y_lo_bc==SOFTWALL_BC and p.pos(1) < plo[1])
            {
                p.pos(1) = two*plo[1] - p.pos(1);
                softwall_arr[i]=1;
            }
            if (y_hi_bc==SOFTWALL_BC and p.pos(1) > phi[1])
            {
                p.pos(1) = two*phi[1] - p.pos(1);
                softwall_arr[i]=1;
            }
            if (z_lo_bc==SOFTWALL_BC and p.pos(2) < plo[2])
            {
                p.pos(2) = two*plo[2] - p.pos(2);
                softwall_arr[i]=1;
            }
            if (z_hi_bc==SOFTWALL_BC and p.pos(2) > phi[2])
            {
                p.pos(2) = two*phi[2] - p.pos(2);
                softwall_arr[i]=1;
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
        auto& soa   = ptile.GetStructOfArrays();
        const size_t np = aos.numParticles();
        ParticleType* pstruct = aos().dataPtr();

        auto pxp_arr = soa.GetRealData(soa_realData::posx_prvs).data();
        auto pyp_arr = soa.GetRealData(soa_realData::posy_prvs).data();
        auto pzp_arr = soa.GetRealData(soa_realData::posz_prvs).data();
        auto softwall_arr = soa.GetIntData(soa_intData::near_softwall).data();

        // now we move the particles
        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            if(softwall_arr[i]==1)
            {
                p.pos(0)=pxp_arr[i];
                p.pos(1)=pyp_arr[i];
                p.pos(2)=pzp_arr[i];
                softwall_arr[i]=0;
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
        auto& soa   = ptile.GetStructOfArrays();
        const size_t np = aos.numParticles();
        ParticleType* pstruct = aos().dataPtr();

        auto vol_arr = soa.GetRealData(soa_realData::volume).data();
        auto total_bridge_vol_arr = soa.GetRealData(soa_realData::total_bridge_volume).data();

        Vector<int*> bridges_vec(MAXBRIDGES);
        for(int br=0; br<MAXBRIDGES; br++) bridges_vec[br] = soa.GetIntData(soa_intData::first_bridge+br).data(); 

        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            // Real MC = (MC_stdev > 0.0) ? std::min(std::max(amrex::RandomNormal(MC_avg, MC_stdev),0.0),0.9):MC_avg;
            Real MC = MC_avg;
	          Real original_lv = p.rdata(aos_realData::liquid_volume);
            p.rdata(aos_realData::liquid_volume) = (MC > FSP || MC == 0.0) ? (p.rdata(aos_realData::density)*vol_arr[i]/liquid_density)*(MC - FSP)/(1.0 - MC):0;
            p.rdata(aos_realData::mass) = p.rdata(aos_realData::density)*vol_arr[i] * (1.0 + MC/(1.0 - MC));
            p.rdata(aos_realData::density) = p.rdata(aos_realData::mass) / vol_arr[i];

	          // Rescale the total liquid participating in bridges
	          if(original_lv > 0.0) total_bridge_vol_arr[i] *= p.rdata(aos_realData::liquid_volume) / original_lv;

	          // If MC < FSP, remove all existing liquid bridges
	          for(int br=0; br<MAXBRIDGES; br++) bridges_vec[br][i] = -1;
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
            //if(p.idata(intData::phase)==0 && ppos[cdir] > cutoff)
            //{
              return(p.rdata(aos_realData::mass));
            //}
            //else
            //{
            //  return(0.0);
            //}
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
            Real pvel[3] = {p.rdata(aos_realData::xvel), p.rdata(aos_realData::yvel), p.rdata(aos_realData::zvel)};
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

    Vector<int> writeflags_real(aos_realData::count+MAXSPECIES-1,0);
    Vector<int> writeflags_int(soa_intData::count,0);

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

    // for(int i=0;i<m_chemptr->nspecies;i++)
    // {
    //     real_data_names[realData::firstspec+i]=m_chemptr->specnames[i];
    // }

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

    // writeflags_real[realData::radius]=1;
    // writeflags_real[realData::xvel]=1;
    // writeflags_real[realData::yvel]=1;
    // writeflags_real[realData::zvel]=1;
    // writeflags_real[realData::fx]=1;
    // writeflags_real[realData::fy]=1;
    // writeflags_real[realData::fz]=1;
    // writeflags_real[realData::xangvel]=1;
    // writeflags_real[realData::yangvel]=1;
    // writeflags_real[realData::zangvel]=1;
    // writeflags_real[realData::taux]=1;
    // writeflags_real[realData::tauy]=1;
    // writeflags_real[realData::tauz]=1;
    // writeflags_real[realData::mass]=1;
    // writeflags_real[realData::temperature]=1;
    // writeflags_real[realData::liquid_volume]=1;
    // writeflags_real[realData::total_bridge_volume]=1;
    // writeflags_real[realData::theta_x]=1;

    // for(int i=0;i<m_chemptr->nspecies;i++)
    // {
    //     writeflags_real[realData::firstspec+i]=1;
    // }
    // if(bonded_sphere_particles){
    //     writeflags_int[intData::type_id] = 1;
    // }

    // WritePlotFile(pltfile, "particles",writeflags_real, 
    //         writeflags_int, real_data_names, int_data_names);

    WritePlotFile(pltfile, "particles");

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

                    if (ls_value < p.rdata(aos_realData::radius)*0.95)
                    {
                        p.id() = -1;
                    }
                });
            }
        }
    }
    Redistribute();
}



