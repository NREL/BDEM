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
        return( std::pow(DEM::k_n/p.rdata(realData::mass)
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
            const int ls_refinement,bool stl_geom_present)
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
    
    //std::map<PairIndex, bool> tile_has_walls;

    /*if(resolve_wall_collisions)
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
            if ((eb_factory != NULL)
                    and ((*flags)[mfi].getType(amrex::grow(bx,1)) == FabType::singlevalued))
            {
                has_wall = true;
            }
            else
            {
                int int_has_wall = 0;
                Real tol = std::min(dx[0], std::min(dx[1], dx[2])) / 2;
                Array4<const Real> const& phi = lsmfab->array(mfi);

#ifdef AMREX_USE_GPU
                Gpu::DeviceScalar<int> has_wall_gpu(int_has_wall);
                int* p_has_wall = has_wall_gpu.dataPtr();
#endif
                amrex::ParallelFor(phibx, [phi,tol,
#ifdef AMREX_USE_GPU
                        p_has_wall]
#else
                        &int_has_wall]
#endif
                        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                        {
                        if(phi(i,j,k) <= tol)
#ifdef AMREX_USE_GPU
                        *p_has_wall = 1;
#else
                        int_has_wall = 1;
#endif
                        });

#ifdef AMREX_USE_GPU
                Gpu::synchronize();
                has_wall = has_wall_gpu.dataValue();
#endif
                has_wall = (int_has_wall > 0);
            }

            tile_has_walls[index] = has_wall;
        }
    }*/

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
            //if (tile_has_walls[index])
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


void BDEMParticleContainer::moveParticles(const amrex::Real& dt,RealVect &gravity,
        int do_chemistry,Real minradfrac)
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
            Real pos_old[3];
            Real pos_new[3];

            pos_old[0]=p.pos(0);
            pos_old[1]=p.pos(1);
            pos_old[2]=p.pos(2);

            p.rdata(realData::xvel) += (p.rdata(realData::fx)/p.rdata(realData::mass) + gravity[XDIR]) * dt;
            p.rdata(realData::yvel) += (p.rdata(realData::fy)/p.rdata(realData::mass) + gravity[YDIR]) * dt;
            p.rdata(realData::zvel) += (p.rdata(realData::fz)/p.rdata(realData::mass) + gravity[ZDIR]) * dt;

            p.pos(0) += p.rdata(realData::xvel) * dt;
            p.pos(1) += p.rdata(realData::yvel) * dt;
            p.pos(2) += p.rdata(realData::zvel) * dt;

            p.rdata(realData::xangvel) += p.rdata(realData::taux) * p.rdata(realData::Iinv) *dt;
            p.rdata(realData::yangvel) += p.rdata(realData::tauy) * p.rdata(realData::Iinv) *dt;
            p.rdata(realData::zangvel) += p.rdata(realData::tauz) * p.rdata(realData::Iinv) *dt;

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

            if (x_lo_bc==HARDWALL_BC and p.pos(0) < plo[0])
            {
                p.pos(0) = two*plo[0] - p.pos(0);
                p.rdata(realData::xvel) = -p.rdata(realData::xvel);
            }
            if (x_hi_bc==HARDWALL_BC and p.pos(0) > phi[0])
            {
                p.pos(0) = two*phi[0] - p.pos(0);
                p.rdata(realData::xvel) = -p.rdata(realData::xvel);
            }
            if (y_lo_bc==HARDWALL_BC and p.pos(1) < plo[1])
            {
                p.pos(1) = two*plo[1] - p.pos(1);
                p.rdata(realData::yvel) = -p.rdata(realData::yvel);
            }
            if (y_hi_bc==HARDWALL_BC and p.pos(1) > phi[1])
            {
                p.pos(1) = two*phi[1] - p.pos(1);
                p.rdata(realData::yvel) = -p.rdata(realData::yvel);
            }
            if (z_lo_bc==HARDWALL_BC and p.pos(2) < plo[2])
            {
                p.pos(2) = two*plo[2] - p.pos(2);
                p.rdata(realData::zvel) = -p.rdata(realData::zvel);
            }
            if (z_hi_bc==HARDWALL_BC and p.pos(2) > phi[2])
            {
                p.pos(2) = two*phi[2] - p.pos(2);
                p.rdata(realData::zvel) = -p.rdata(realData::zvel);
            }

        });

    }
}

void BDEMParticleContainer::writeParticles(const int n)
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
    real_data_names.push_back("volume");
    real_data_names.push_back("mass");
    real_data_names.push_back("density");
    real_data_names.push_back("temperature");

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
    writeflags_real[realData::temperature]=1;
    for(int i=0;i<m_chemptr->nspecies;i++)
    {
        writeflags_real[realData::firstspec+i]=1;
    }

    WritePlotFile(pltfile, "particles",writeflags_real, 
            writeflags_int, real_data_names, int_data_names);

}
