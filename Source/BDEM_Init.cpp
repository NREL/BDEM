#include <BDEM_ParticleContainer.H>
#include <BDEM_Collide_Utils.H>
#include <stl_tools/STLtools.H>

void BDEMParticleContainer::InitParticles (const std::string& filename,bool &do_heat_transfer)
{

    // only read the file on the IO proc
    if (ParallelDescriptor::IOProcessor())  
    {
        std::ifstream ifs;
        ifs.open(filename.c_str(), std::ios::in);

        if (!ifs.good())
        {
            amrex::FileOpenFailed(filename);
        }

        int np = -1;
        ifs >> np >> std::ws;

        // Issue an error if nparticles = 0 is specified
        if ( np == -1 )
        {
            Abort("\nCannot read number of particles from particle_input.dat: file is corrupt.\
                    \nPerhaps you forgot to specify the number of particles on the first line??? ");
        }

        // we add all the particles to grid 0 and tile 0 and let
        // Redistribute() put them in the right places.
        const int lev  = 0;
        const int grid = 0;
        const int tile = 0;

        Gpu::HostVector<ParticleType> host_particles;
        
        auto& particle_tile = DefineAndReturnParticleTile(lev,grid,tile);

        for (int i = 0; i < np; i++) 
        {
            ParticleType p;
            // Set id and cpu for this particle
            p.id()  = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

            // Read from input file
            ifs >> p.idata(intData::phase);
            ifs >> p.pos(0);
            ifs >> p.pos(1);
            ifs >> p.pos(2);
            ifs >> p.rdata(realData::radius);
            ifs >> p.rdata(realData::density);
            ifs >> p.rdata(realData::E);
            ifs >> p.rdata(realData::nu);
            ifs >> p.rdata(realData::xvel);
            ifs >> p.rdata(realData::yvel);
            ifs >> p.rdata(realData::zvel);

            if(do_heat_transfer)
            {
                ifs >> p.rdata(realData::temperature);
            }
            else
            {
               p.rdata(realData::temperature)=NTP_TEMP;
            }

            //set initial radius
            p.rdata(realData::radinit)=p.rdata(realData::radius);

            // Set other particle properties
            p.rdata(realData::volume)      = fourbythree*PI*pow(p.rdata(realData::radius),three);
            p.rdata(realData::mass)        = p.rdata(realData::density)*p.rdata(realData::volume);
            p.rdata(realData::Iinv)        = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
            p.rdata(realData::xangvel)     = amrex::Random();
            p.rdata(realData::yangvel)     = amrex::Random();
            p.rdata(realData::zangvel)     = amrex::Random();

            p.rdata(realData::fx) = zero;
            p.rdata(realData::fy) = zero;
            p.rdata(realData::fz) = zero;
            p.rdata(realData::taux) = zero;
            p.rdata(realData::tauy) = zero;
            p.rdata(realData::tauz) = zero;

            
            //FIXME: get chemistry data from inputs file
            for(int sp=0;sp<MAXSPECIES;sp++)
            {
                p.rdata(realData::firstspec+sp)=0.0;
            }

            // Add everything to the data structure
            host_particles.push_back(p);

            if (!ifs.good())
            {
                amrex::Abort("Error initializing particles from Ascii file. \n");
            }
        }
        
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);
    }
    
    Redistribute();
}

void BDEMParticleContainer::InitChemSpecies(amrex::Real Yis[MAXSPECIES])
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    const auto dx = Geom(lev).CellSizeArray();
    auto& plev  = GetParticles(lev);

    amrex::Real Yi_captured[MAXSPECIES];

    for(int i=0;i<MAXSPECIES;i++)
    {
        Yi_captured[i]=Yis[i];
    }

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

                for(int i=0;i<MAXSPECIES;i++)
                {
                    //p.rdata(realData::firstspec+i)=Yis[i];
                    p.rdata(realData::firstspec+i)=Yi_captured[i];
                }
        });
    }
}

void BDEMParticleContainer::InitChemSpecies(int ndomains, Real *mincoords,
        Real *maxcoords,Real *spec_massfracs)
{
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    const auto dx = Geom(lev).CellSizeArray();
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

            Real x=p.pos(0);
            Real y=p.pos(1);
            Real z=p.pos(2);

            for(int d=0;d<ndomains;d++)
            {
                if(x>mincoords[d*AMREX_SPACEDIM+0] && x<=maxcoords[d*AMREX_SPACEDIM+0] &&
                        y>mincoords[d*AMREX_SPACEDIM+1] && y<=maxcoords[d*AMREX_SPACEDIM+1] &&
                        z>mincoords[d*AMREX_SPACEDIM+2] && z<=maxcoords[d*AMREX_SPACEDIM+2])
                {
                    for(int i=0;i<MAXSPECIES;i++)
                    {
                        p.rdata(realData::firstspec+i)=spec_massfracs[d*MAXSPECIES+i];
                    }

                }

            }

        });
    }
}

void BDEMParticleContainer::removeParticlesOutsideBoundary(const MultiFab *lsmfab,
        const EBFArrayBoxFactory *ebfactory,
        const int ls_refinement)
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);
    const FabArray<EBCellFlagFab>* flags = &(ebfactory->getMultiEBCellFlagFab());
    const Geometry& geom = Geom(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

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

        if ((*flags)[mfi].getType(bx) == FabType::covered)
        {
            amrex::ParallelFor(np,[=]
            AMREX_GPU_DEVICE (int i) noexcept
            {
                ParticleType& p = pstruct[i];
                p.id()=-1;   
            });
        }
        else
        {
            const auto phiarr = lsmfab->array(mfi);
            amrex::ParallelFor(np,[=]
                  AMREX_GPU_DEVICE (int i) noexcept
            {
                ParticleType& p = pstruct[i];
                Real rp = p.rdata(realData::radius);

                Real ls_value = get_levelset_value(p, ls_refinement, phiarr, plo, dx);
                if(ls_value < 0.0)
                {
                    p.id()=-1;   
                }
            });
        }
    }
    Redistribute();
}

void BDEMParticleContainer::removeParticlesInsideSTL(Vector<Real> outside_point)
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);
    const Geometry& geom = Geom(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    Array<Real,AMREX_SPACEDIM> po_arr{outside_point[0],outside_point[1],outside_point[2]};

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
            Real ploc[3]={p.pos(0),p.pos(1),p.pos(2)};
            Real ploc_t[3]={0.0};
            Real t1[3],t2[3],t3[3];
            Real outp[]={po_arr[0],po_arr[1],po_arr[2]};
            int num_intersects=0;

            for(int dim=0;dim<3;dim++)
            {
                for(int j=0;j<3;j++)
                {
                    ploc_t[dim] += STLtools::eigdirs[3*dim+j]*ploc[j];
                }
            }

            if( (ploc_t[0]>STLtools::bbox_lo[0]) && 
               (ploc_t[0]<STLtools::bbox_hi[0]) &&
               (ploc_t[1]>STLtools::bbox_lo[1]) &&
               (ploc_t[1]<STLtools::bbox_hi[1]) &&
               (ploc_t[2]>STLtools::bbox_lo[2]) &&
               (ploc_t[2]<STLtools::bbox_hi[2]) )
            {
                for(int tr=0;tr<STLtools::num_tri;tr++)
                {
                    t1[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+0];
                    t1[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+1];
                    t1[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+2];

                    t2[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+3];
                    t2[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+4];
                    t2[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+5];

                    t3[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+6];
                    t3[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+7];
                    t3[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+8];

                    num_intersects += (1-STLtools::lineseg_tri_intersect(outp,ploc,t1,t2,t3));
                }
                if(num_intersects%2 == 1)
                {
                    p.id()=-1;   
                }
            }
        });
    }

    Redistribute();
}

void BDEMParticleContainer::checkParticlesInsideSTL(Vector<Real> outside_point)
{
    const int lev = 0;
    auto& plev  = GetParticles(lev);
    const Geometry& geom = Geom(lev);
    const auto dxi = geom.InvCellSizeArray();
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    Array<Real,AMREX_SPACEDIM> po_arr{outside_point[0],outside_point[1],outside_point[2]};

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
            Real ploc[3]={p.pos(0),p.pos(1),p.pos(2)};
            Real ploc_t[3]={0.0};
            Real t1[3],t2[3],t3[3];
            Real outp[]={po_arr[0],po_arr[1],po_arr[2]};
            int num_intersects=0;

            for(int dim=0;dim<3;dim++)
            {
                for(int j=0;j<3;j++)
                {
                    ploc_t[dim] += STLtools::eigdirs[3*dim+j]*ploc[j];
                }
            }

            if( (ploc_t[0]>STLtools::bbox_lo[0]) && 
               (ploc_t[0]<STLtools::bbox_hi[0]) &&
               (ploc_t[1]>STLtools::bbox_lo[1]) &&
               (ploc_t[1]<STLtools::bbox_hi[1]) &&
               (ploc_t[2]>STLtools::bbox_lo[2]) &&
               (ploc_t[2]<STLtools::bbox_hi[2]) )
            {
                for(int tr=0;tr<STLtools::num_tri;tr++)
                {
                    t1[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+0];
                    t1[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+1];
                    t1[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+2];

                    t2[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+3];
                    t2[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+4];
                    t2[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+5];

                    t3[0]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+6];
                    t3[1]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+7];
                    t3[2]=STLtools::tri_pts[tr*STLtools::ndata_per_tri+8];

                    num_intersects += (1-STLtools::lineseg_tri_intersect(outp,ploc,t1,t2,t3));
                }
                if(num_intersects%2 == 1)
                {
#ifndef AMREX_USE_GPU
                    amrex::Print()<<"particle inside stl:"<<ploc[0]<<"\t"<<ploc[1]<<"\t"<<ploc[2]<<"\n";
#else
                    amrex::Abort("particle inside stl");
#endif

                }
            }
        });
    }

}


void BDEMParticleContainer::InitParticles (Real mincoords[THREEDIM],Real maxcoords[THREEDIM], 
                                           Real meanvel[THREEDIM], Real fluctuation[THREEDIM], Real rad, Real dens,
                                           Real E, Real nu, Real temp,
                                           Real spec[MAXSPECIES],
                                           int do_multi_part_per_cell)
{
    int lev = 0;
    Real x,y,z,x0,y0,z0;

    Real dx = Geom(lev).CellSize(0);
    Real dy = Geom(lev).CellSize(1);
    Real dz = Geom(lev).CellSize(2);
    const Real* plo = Geom(lev).ProbLo();

    std::mt19937 mt(0451);
    std::uniform_real_distribution<double> dist(0.4, 0.6);


    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) 
    {

        const Box& tile_box = mfi.tilebox();
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particle_tile = GetParticles(lev)[std::make_pair(grid_id, tile_id)];

        Gpu::HostVector<ParticleType> host_particles;

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) 
        {
            if(do_multi_part_per_cell == 0)
            {
                x = plo[XDIR] + (iv[XDIR] + dist(mt))*dx;
                y = plo[YDIR] + (iv[YDIR] + dist(mt))*dy;
                z = plo[ZDIR] + (iv[ZDIR] + dist(mt))*dz;

                if(x>=mincoords[XDIR] && x<=maxcoords[XDIR] &&
                   y>=mincoords[YDIR] && y<=maxcoords[YDIR] &&
                   z>=mincoords[ZDIR] && z<=maxcoords[ZDIR])
                {
                    ParticleType p = generate_particle(x,y,z,
                                                       meanvel[XDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                       meanvel[YDIR] + fluctuation[YDIR]*(amrex::Random()-half),
                                                       meanvel[ZDIR] + fluctuation[ZDIR]*(amrex::Random()-half),
                                                       dens, rad, E, nu, 
                                                       temp, spec);
                    host_particles.push_back(p);
                }
            }
            else
            {
                x0 = plo[XDIR]+iv[XDIR]*dx;
                y0 = plo[YDIR]+iv[YDIR]*dy;
                z0 = plo[ZDIR]+iv[ZDIR]*dz;

                for(int k=0;k<2;k++)
                {
                    for(int j=0;j<2;j++)
                    {
                        for(int i=0;i<2;i++)
                        {
                            //x = x0 + (i+dist(mt))*half*dx;
                            //y = y0 + (j+dist(mt))*half*dy;
                            //z = z0 + (k+dist(mt))*half*dz;
                            x = x0 + (i+half)*half*dx;
                            y = y0 + (j+half)*half*dy;
                            z = z0 + (k+half)*half*dz;

                            if(x>=mincoords[XDIR] and x<=maxcoords[XDIR] and 
                               y>=mincoords[YDIR] and y<=maxcoords[YDIR] and
                               z>=mincoords[ZDIR] and z<=maxcoords[ZDIR])
                            {
                                ParticleType p = generate_particle(x,y,z,
                                                                   meanvel[XDIR] + fluctuation[XDIR]*(amrex::Random()-half),
                                                                   meanvel[YDIR] + fluctuation[YDIR]*(amrex::Random()-half),
                                                                   meanvel[ZDIR] + fluctuation[ZDIR]*(amrex::Random()-half),
                                                                   dens, rad, E, nu, 
                                                                   temp, spec);
                                host_particles.push_back(p);
                            }
                        }
                    } 
                }
            }
        }

        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(Gpu::hostToDevice,
                  host_particles.begin(),
                  host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);

    }

    // We shouldn't need this if the particles are tiled with one tile per grid, but otherwise
    // we do need this to move particles from tile 0 to the correct tile.
    Redistribute();
}

BDEMParticleContainer::ParticleType BDEMParticleContainer::generate_particle(Real x,Real y,Real z,
                                                                             Real velx, Real vely, Real velz,
                                                                             Real dens, Real rad, Real E, Real nu,
                                                                             Real temp, Real spec[MAXSPECIES])
{
    ParticleType p;
    p.id()  = ParticleType::NextID();
    p.cpu() = ParallelDescriptor::MyProc();                

    p.pos(XDIR) = x;
    p.pos(YDIR) = y;
    p.pos(ZDIR) = z;

    p.idata(intData::phase) = 0;
    p.rdata(realData::radius) = rad;
    p.rdata(realData::radinit) = rad;
    p.rdata(realData::E) = E;
    p.rdata(realData::nu) = nu;

    p.rdata(realData::density) = dens;
    p.rdata(realData::xvel) = velx;
    p.rdata(realData::yvel) = vely;
    p.rdata(realData::zvel) = velz;
    p.rdata(realData::temperature)=temp;

    // Set other particle properties
    p.rdata(realData::volume)      = fourbythree*PI*pow(p.rdata(realData::radius),three);
    p.rdata(realData::mass)        = p.rdata(realData::density)*p.rdata(realData::volume);
    p.rdata(realData::Iinv)        = 2.5/(p.rdata(realData::mass)*pow(p.rdata(realData::radius),two));
    p.rdata(realData::xangvel)     = zero;
    p.rdata(realData::yangvel)     = zero;
    p.rdata(realData::zangvel)     = zero;

    p.rdata(realData::fx) = zero;
    p.rdata(realData::fy) = zero;
    p.rdata(realData::fz) = zero;
    p.rdata(realData::taux) = zero;
    p.rdata(realData::tauy) = zero;
    p.rdata(realData::tauz) = zero;

    for(int i=0;i<MAXSPECIES;i++)
    {
        p.rdata(realData::firstspec+i)=spec[i];
    }

    return(p);
}
