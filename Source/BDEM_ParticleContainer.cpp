#include <BDEM_ParticleContainer.H>
#include <BDEM_Collide_Utils.H>
#include <stl_tools/STLtools.H>
#include <MoveUtils.H>
#include <WallTemp.H>

using namespace amrex;
using namespace std;

Real BDEMParticleContainer::compute_coll_timescale()
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
        return( std::pow(DEM::k_n/(p.rdata(realData::mass)/p.idata(intData::num_comp_sphere))
                *PI*PI/(PI*PI + log(DEM::e_n)*log(DEM::e_n))
                ,-half)
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
            const int glued_sphere_particles,
            const int bonded_sphere_particles,
            const int liquid_bridging, 
            const Real init_force, const int init_force_dir, const int init_force_comp,
            const Real cb_force, const int cb_dir)
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


void BDEMParticleContainer::moveParticles(const amrex::Real& dt,
        int do_chemistry,Real minradfrac, const int glued_sphere_particles)
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

              p.rdata(realData::xvel) += (p.rdata(realData::fx)/p.rdata(realData::mass)) * dt;
              p.rdata(realData::yvel) += (p.rdata(realData::fy)/p.rdata(realData::mass)) * dt;
              p.rdata(realData::zvel) += (p.rdata(realData::fz)/p.rdata(realData::mass)) * dt;

              p.pos(0) += p.rdata(realData::xvel) * dt;
              p.pos(1) += p.rdata(realData::yvel) * dt;
              p.pos(2) += p.rdata(realData::zvel) * dt;

              if(glued_sphere_particles){
                  // Angular velocity is updated in body-fixed frame of reference so that diagonal inertia tensor can be used
                  // I d (w_body)/dt = -w_body x (w_body I) + R tau
  
                  // Calculate the body-fixed angular velocity
                  Real angvel_inert[THREEDIM] = {p.rdata(realData::xangvel), p.rdata(realData::yangvel), p.rdata(realData::zangvel)};
                  Real angvel_body[THREEDIM];
                  rotate_vector_to_body(p, angvel_inert, angvel_body);

                  // Calculate the cross product term
                  Real wI[THREEDIM];
                  Real cpdt[THREEDIM];
                  wI[XDIR] = angvel_body[XDIR] / p.rdata(realData::Ixinv);
                  wI[YDIR] = angvel_body[YDIR] / p.rdata(realData::Iyinv);
                  wI[ZDIR] = angvel_body[ZDIR] / p.rdata(realData::Izinv);
                  crosspdt(angvel_body, wI, cpdt);

                  // Rotate the torque vector to body-fixed frame
                  Real tau_inert[THREEDIM] = {p.rdata(realData::taux), p.rdata(realData::tauy), p.rdata(realData::tauz)};
                  Real tau_body[THREEDIM];
                  rotate_vector_to_body(p, tau_inert, tau_body);

                  // Update the angular velocity in the body-fixed reference frame
                  angvel_body[XDIR] += (tau_body[XDIR] - cpdt[XDIR]) * p.rdata(realData::Ixinv) * dt;
                  angvel_body[YDIR] += (tau_body[YDIR] - cpdt[YDIR]) * p.rdata(realData::Iyinv) * dt;
                  angvel_body[ZDIR] += (tau_body[ZDIR] - cpdt[ZDIR]) * p.rdata(realData::Izinv) * dt;

                  // Update the quaternion components with the updated angular velocity
                  Real q0 = p.rdata(realData::q0);
                  Real q1 = p.rdata(realData::q1);
                  Real q2 = p.rdata(realData::q2);
                  Real q3 = p.rdata(realData::q3);

                  Real dq0 = (dt/2.0) * (-q1*angvel_inert[XDIR] - q2*angvel_inert[YDIR] - q3*angvel_inert[ZDIR]);
                  Real dq1 = (dt/2.0) * ( q0*angvel_inert[XDIR] + q3*angvel_inert[YDIR] - q2*angvel_inert[ZDIR]);
                  Real dq2 = (dt/2.0) * (-q3*angvel_inert[XDIR] + q0*angvel_inert[YDIR] + q1*angvel_inert[ZDIR]);
                  Real dq3 = (dt/2.0) * ( q2*angvel_inert[XDIR] - q1*angvel_inert[YDIR] + q0*angvel_inert[ZDIR]);

                  q0 += dq0;
                  q1 += dq1;
                  q2 += dq2;
                  q3 += dq3;

                  // Normalize quaternion components to ensure sum_i (q_i)^2 = 1
                  Real qmag = sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
                  
                  if(qmag > TINYVAL){
                      p.rdata(realData::q0) = (q0 - amrex::Math::copysign(1.0,q0)*TINYVAL)/qmag;
                      p.rdata(realData::q1) = (q1 - amrex::Math::copysign(1.0,q1)*TINYVAL)/qmag;
                      p.rdata(realData::q2) = (q2 - amrex::Math::copysign(1.0,q2)*TINYVAL)/qmag;
                      p.rdata(realData::q3) = (q3 - amrex::Math::copysign(1.0,q3)*TINYVAL)/qmag;
                  }

                  // Make sure single particles are not rotated
                  if(p.idata(intData::num_comp_sphere) == 1){
                      p.rdata(realData::q0) = 1.0;
                      p.rdata(realData::q1) = zero;
                      p.rdata(realData::q2) = zero;
                      p.rdata(realData::q3) = zero;
                  }

                  // Rotate updated body-fixed angular velocity back to inertial frame
                  rotate_vector_to_inertial(p, angvel_body, angvel_inert);
                  p.rdata(realData::xangvel) = angvel_inert[XDIR];
                  p.rdata(realData::yangvel) = angvel_inert[YDIR];
                  p.rdata(realData::zangvel) = angvel_inert[ZDIR];

                  // Calculate principal axis components in inertial reference frame
                  Real pa_body[THREEDIM] = {1.0, 0.0, 0.0};
                  Real pa_inert[THREEDIM];
                  rotate_vector_to_inertial(p, pa_body, pa_inert);
                  p.rdata(realData::pax) = pa_inert[XDIR];
                  p.rdata(realData::pay) = pa_inert[YDIR];
                  p.rdata(realData::paz) = pa_inert[ZDIR];
              } else{
                  p.rdata(realData::xangvel) += p.rdata(realData::taux) * p.rdata(realData::Iinv) *dt;
                  p.rdata(realData::yangvel) += p.rdata(realData::tauy) * p.rdata(realData::Iinv) *dt;
                  p.rdata(realData::zangvel) += p.rdata(realData::tauz) * p.rdata(realData::Iinv) *dt;
              }

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

void BDEMParticleContainer::computeMoistureContent(Real moisture_content, Real contact_angle, Real liquid_density, Real fiber_sat_pt)
{
    BL_PROFILE("BDEMParticleContainer::computeMoistureContent");

    Real MC = moisture_content/100.0;
    Real FSP = fiber_sat_pt/100.0;

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

        // now we move the particles
        amrex::ParallelFor(np,[=]
        AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p = pstruct[i];
            p.rdata(realData::liquid_volume) = (MC > FSP) ? (p.rdata(realData::density)*p.rdata(realData::volume)/liquid_density)*(MC - FSP)/(1 - MC):0;
            p.rdata(realData::mass) = p.rdata(realData::density)*p.rdata(realData::volume) * (1.0 + MC/(1.0 - MC));
        });
    }
}

void BDEMParticleContainer::writeParticles(const int n, const int glued_sphere_particles)
{
    BL_PROFILE("BDEMParticleContainer::writeParticles");
    const std::string& pltfile = amrex::Concatenate("plt", n, 5);

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
    real_data_names.push_back("Ixinv");
    real_data_names.push_back("Iyinv");
    real_data_names.push_back("Izinv");
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
    real_data_names.push_back("q0");
    real_data_names.push_back("q1");
    real_data_names.push_back("q2");
    real_data_names.push_back("q3");
    real_data_names.push_back("pax");
    real_data_names.push_back("pay");
    real_data_names.push_back("paz");
    real_data_names.push_back("liquid_volume");
    real_data_names.push_back("total_bridge_volume");
    real_data_names.push_back("fx_bond");
    real_data_names.push_back("fy_bond");
    real_data_names.push_back("fz_bond");
    real_data_names.push_back("taux_bond");
    real_data_names.push_back("tauy_bond");
    real_data_names.push_back("tauz_bond");

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
    int_data_names.push_back("num_comp_sphere");

    for(int i=0;i<MAXBRIDGES;i++)
    {
        std::string bridgeidx = amrex::Concatenate("p2idx",i,2);
        int_data_names.push_back(bridgeidx);
        bridgeidx = amrex::Concatenate("p1cs",i,2);
        int_data_names.push_back(bridgeidx);
        bridgeidx = amrex::Concatenate("p2cs",i,2);
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
    writeflags_real[realData::Ixinv]=1;
    writeflags_real[realData::Iyinv]=1;
    writeflags_real[realData::Izinv]=1;
    writeflags_real[realData::q0]=1;
    writeflags_real[realData::q1]=1;
    writeflags_real[realData::q2]=1;
    writeflags_real[realData::q3]=1;
    writeflags_real[realData::temperature]=1;
    writeflags_real[realData::liquid_volume]=1;
    writeflags_real[realData::total_bridge_volume]=1;

    // for debug
    // writeflags_real[realData::overlap_n]=1;
    // writeflags_real[realData::kn]=1;
    // writeflags_real[realData::kt]=1;
    // writeflags_real[realData::eta_n]=1;
    // writeflags_real[realData::eta_t]=1;
    // writeflags_real[realData::overlap_n_y]=1;
    // writeflags_real[realData::fn_y]=1;
    // writeflags_real[realData::vr_n_y]=1;
    // writeflags_real[realData::overlap_t_y]=1;
    // writeflags_real[realData::ft_y]=1;
    // writeflags_real[realData::vr_t_y]=1;

    for(int i=0;i<m_chemptr->nspecies;i++)
    {
        writeflags_real[realData::firstspec+i]=1;
    }
    if(glued_sphere_particles){
        writeflags_real[realData::pax] = 1;
        writeflags_real[realData::pay] = 1;
        writeflags_real[realData::paz] = 1;
        writeflags_int[intData::num_comp_sphere] = 1;
    }

    WritePlotFile(pltfile, "particles",writeflags_real, 
            writeflags_int, real_data_names, int_data_names);

}

void BDEMParticleContainer::createGluedSpheres(BDEMParticleContainer& pin)
{
    // Extract particle tile from input BPC
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev_in = pin.GetParticles(lev);
    auto& plev_out = GetParticles(lev);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile_in = plev_in[index];
        auto& ptile_out = plev_out[index];
        auto& aos_in   = ptile_in.GetArrayOfStructs();
        const size_t np = aos_in.numParticles();

        ParticleType* pstruct_in = aos_in().dataPtr();
        Gpu::HostVector<ParticleType> host_particles;

        for(int i = 0; i<np; i++){
            ParticleType& p_in = pstruct_in[i];
            for(int pc = 0; pc < p_in.idata(intData::num_comp_sphere); pc++){
                ParticleType pcomp;

                // Calculate inertial frame position of each component sphere
                Real ppos_inert[THREEDIM];
                Real ppos_body[THREEDIM] = {p_in.rdata(realData::radius)*((p_in.idata(intData::num_comp_sphere) - 1) - 2.0*pc), 0.0, 0.0};
                rotate_vector_to_inertial(p_in, ppos_body, ppos_inert);
                Real pposx = p_in.pos(0) + ppos_inert[XDIR];
                Real pposy = p_in.pos(1) + ppos_inert[YDIR];
                Real pposz = p_in.pos(2) + ppos_inert[ZDIR];

                pcomp.id()   = ParticleType::NextID();
                pcomp.cpu()  = ParallelDescriptor::MyProc();
                
                pcomp.pos(0) = pposx;
                pcomp.pos(1) = pposy;
                pcomp.pos(2) = pposz;

                pcomp.rdata(realData::radius) = p_in.rdata(realData::radius);
                pcomp.rdata(realData::radinit) = p_in.rdata(realData::radinit);
                pcomp.rdata(realData::density) = p_in.rdata(realData::density);
                pcomp.rdata(realData::temperature) = p_in.rdata(realData::temperature);
                pcomp.rdata(realData::mass) = p_in.rdata(realData::mass) / p_in.idata(intData::num_comp_sphere);
                pcomp.rdata(realData::volume) = p_in.rdata(realData::volume) / p_in.idata(intData::num_comp_sphere);

                Real rvel[THREEDIM];
                Real avel[THREEDIM] = {p_in.rdata(realData::xangvel), p_in.rdata(realData::yangvel), p_in.rdata(realData::zangvel)};
                crosspdt(avel, ppos_inert, rvel);
                pcomp.rdata(realData::xvel) = p_in.rdata(realData::xvel) + rvel[XDIR];
                pcomp.rdata(realData::yvel) = p_in.rdata(realData::yvel) + rvel[YDIR];
                pcomp.rdata(realData::zvel) = p_in.rdata(realData::zvel) + rvel[ZDIR];

                pcomp.idata(intData::num_comp_sphere) = p_in.idata(intData::num_comp_sphere);
                pcomp.rdata(realData::euler_angle_x) = p_in.rdata(realData::euler_angle_x);
                pcomp.rdata(realData::euler_angle_y) = p_in.rdata(realData::euler_angle_y);
                pcomp.rdata(realData::euler_angle_z) = p_in.rdata(realData::euler_angle_z);
                pcomp.rdata(realData::q0) = p_in.rdata(realData::q0);
                pcomp.rdata(realData::q1) = p_in.rdata(realData::q1);
                pcomp.rdata(realData::q2) = p_in.rdata(realData::q2);
                pcomp.rdata(realData::q3) = p_in.rdata(realData::q3);
                pcomp.rdata(realData::pax) = p_in.rdata(realData::pax);
                pcomp.rdata(realData::pay) = p_in.rdata(realData::pay);
                pcomp.rdata(realData::paz) = p_in.rdata(realData::paz);
                pcomp.rdata(realData::Iinv) = p_in.rdata(realData::Iinv);
                pcomp.rdata(realData::Ixinv) = p_in.rdata(realData::Ixinv);
                pcomp.rdata(realData::Iyinv) = p_in.rdata(realData::Iyinv);
                pcomp.rdata(realData::Izinv) = p_in.rdata(realData::Izinv);
                pcomp.rdata(realData::xangvel) = p_in.rdata(realData::xangvel);
                pcomp.rdata(realData::yangvel) = p_in.rdata(realData::yangvel);
                pcomp.rdata(realData::zangvel) = p_in.rdata(realData::zangvel);
                pcomp.rdata(realData::fx) = p_in.rdata(realData::fx);
                pcomp.rdata(realData::fy) = p_in.rdata(realData::fy);
                pcomp.rdata(realData::fz) = p_in.rdata(realData::fz);
                pcomp.rdata(realData::taux) = p_in.rdata(realData::taux);
                pcomp.rdata(realData::tauy) = p_in.rdata(realData::tauy);
                pcomp.rdata(realData::tauz) = p_in.rdata(realData::tauz);
                pcomp.rdata(realData::liquid_volume) = p_in.rdata(realData::liquid_volume) / p_in.idata(intData::num_comp_sphere);
                pcomp.rdata(realData::total_bridge_volume) = p_in.rdata(realData::total_bridge_volume) / p_in.idata(intData::num_comp_sphere);
                pcomp.rdata(realData::fx_bond) = p_in.rdata(realData::fx_bond);
                pcomp.rdata(realData::fy_bond) = p_in.rdata(realData::fy_bond);
                pcomp.rdata(realData::fz_bond) = p_in.rdata(realData::fz_bond);
                pcomp.rdata(realData::taux_bond) = p_in.rdata(realData::taux_bond);
                pcomp.rdata(realData::tauy_bond) = p_in.rdata(realData::tauy_bond);
                pcomp.rdata(realData::tauz_bond) = p_in.rdata(realData::tauz_bond);

                for(int br=0; br<MAXBRIDGES; br++){
                    pcomp.idata(intData::first_bridge+3*br) = -1;
                    pcomp.idata(intData::first_bridge+3*br+1) = -1;
                    pcomp.idata(intData::first_bridge+3*br+2) = -1;
                }
                for(int sp=0;sp<MAXSPECIES;sp++)
                {
                    pcomp.rdata(realData::firstspec+sp)=p_in.rdata(realData::firstspec+sp);
                }
                for(int b=0; b<MAXBONDS; b++) pcomp.idata(intData::first_bond + b) = -1;
   
                host_particles.push_back(pcomp);
            }
        }

        auto old_size = ptile_out.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        ptile_out.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  ptile_out.GetArrayOfStructs().begin() + old_size);
    }
}



